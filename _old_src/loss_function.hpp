#ifndef LOSS_FUNCTION_HPP_
#define LOSS_FUNCTION_HPP_

#include <vector>
#include <unordered_set>
#include <optional>
#include <expected>

#include <assert.h>
#include <omp.h>
#include <stdio.h>

#include "core_convertor.hpp"
#include "utils.hpp"

#include "parallel_simulation.hpp"


struct SimulationResult {
    std::unordered_map<SpeciesId, double> old_values;
    // TODO: deprecated constrained
    std::unordered_map<SpeciesId, double> constrained;
    std::unordered_map<SpeciesId, double> not_constrained;
    std::vector<std::pair<SpeciesId, double>> mins;
    std::vector<std::pair<SpeciesId, double>> maxs;
};

struct Fitness {
    // TODO: padding
    int simulation_id;
    double fitness;
    bool error;

    Fitness(int id, double fitness): simulation_id(id), fitness(fitness), error(false) { }

    Fitness(int id, bool error): simulation_id(id), error(error) {}

    Fitness(int id, bool error, double fitness): simulation_id(id), fitness(fitness), error(error) {}
};

bool comp(const std::pair<SpeciesId, double>& a, const std::pair<SpeciesId, double>& b) {
    return a.second < b.second;
}

class ErrorHandler {
    std::unordered_set<SpeciesId> constrained_species;
    std::vector<std::pair<SpeciesId, double>> real_conc;
    std::unordered_set<SpeciesId> outputs;
protected:

    const std::vector<std::pair<SpeciesId, double>>& get_real_concs() const {
        return this->real_conc;
    }

    const std::unordered_set<SpeciesId>& get_outputs() const {
        return this->outputs;
    }

public:
    virtual ~ErrorHandler() = default;
    virtual double error(std::expected<SimulationResult,double> sim_result) const = 0;
    virtual void add_output(SpeciesId id) {
        this->outputs.insert(id);
    }

    double get_real_conc(SpeciesId id) const {
        for (const auto& pair : this->real_conc) {
            if (pair.first == id) {
                return pair.second;
            }
        }
        throw std::runtime_error("SpeciesId not found in real_conc");
    }

    void add_real_concentration(const char *id, double value) {
        control(constrained_species.insert(id).second);
        this->real_conc.push_back(std::pair(id, value));
    }

    void order_real_concentration() {
        std::sort(this->real_conc.begin(), this->real_conc.end(), comp);
    }

    const std::unordered_set<SpeciesId>& get_constrained_species() const {
        return this->constrained_species;
    }

    double from_trajectory(const std::vector<std::pair<SpeciesId, std::vector<double>>> *trajectory, double x) {
        std::unordered_map<SpeciesId, double> species_avg_values;
        std::unordered_map<SpeciesId, double> old_species_avg_values;
        std::vector<std::pair<SpeciesId, double>> mins;
        std::vector<std::pair<SpeciesId, double>> maxs;

        const std::string prefix = "avg_";
        for (const auto& [species_id, values] : *trajectory) {
            if (species_id.compare(0, prefix.size(), prefix) == 0) {
                if (!values.empty()) {
                    species_avg_values[species_id] = values[values.size() - 1];
                    int eightyPercentRow = static_cast<int>(values.size() * x);
                    if (eightyPercentRow < 0) eightyPercentRow = 0;
                    if (eightyPercentRow >= static_cast<int>(values.size())) eightyPercentRow = static_cast<int>(values.size()) - 1;
                    old_species_avg_values[species_id] = values[eightyPercentRow];
                }
            } else {
                double min = values[0];
                double max = values[0];
                for (double value : values) {
                    if (value < min) {
                        min = value;
                    }
                    if (value > max) {
                        max = value;
                    }
                }
                mins.emplace_back(species_id, min);
                maxs.emplace_back(species_id, max);
            }
        }

        SimulationResult result;

        result.old_values = old_species_avg_values;
        result.not_constrained = species_avg_values;
        result.mins = mins;
        result.maxs = maxs;

        return this->error(std::expected<SimulationResult,double>(std::move(result)));
    }
};

#define WATER_CONC 0.5

class OrderingError: public ErrorHandler {
public:

    ~OrderingError() override = default;

    OrderingError() = default;

    double error(std::expected<SimulationResult,double> sim_result) const override {
        /* double err = 0.0;
        if(sim_result) {
            SimulationResult simResult = *sim_result;
            std::vector<std::pair<SpeciesId,double>>& result = simResult.constrained;
            std::sort(result.begin(), result.end(), comp);
            auto real_conc = this->get_real_concs();
            control(result.size() == real_conc.size());
            for (size_t i = 0; i < result.size(); ++i) {
                bool not_good = false;    
                if(result[i].second > WATER_CONC || result[i].second < -1e-6) {
                    // unstable
                    not_good = true;
                }
                if(result[i].first != real_conc[i].first || not_good) {
                    double a = std::log10(std::abs((result[i].second + LITTLE_EPSILON) / (this->get_real_conc(result[i].first)+ LITTLE_EPSILON)));
                    if(not_good) {
                        err += 1e4*a*a;
                    } else {
                        err += a * a;
                    }
                }
            }

            const auto& outputs = this->get_outputs();
            
            for (const auto& nc : simResult.not_constrained) {
                if (outputs.find(nc.first) == outputs.end() && (nc.second > WATER_CONC || nc.second < -1e-6)) {
                    // not bounded
                    double a = std::log10(std::abs((nc.second + LITTLE_EPSILON)));
                    err += a;
                }
            } 


        } else {
            // TODO: scegli l'orizzonte di simulazione
            double error = sim_result.error();
            if(error < 1e-12) {
                error = 1e-12;
            }
            err = 1e6*(1/((error)/100.0));
        }
        return err; */
        TODO("sistemare");
    }

    
};

class ClassicalError: public ErrorHandler {

    double S = 1.5;
    double p = 0.8;
public:
    ~ClassicalError() override = default;

    ClassicalError() = default;

    ClassicalError(double S, double p): S(S), p(p) {}

    double error(std::expected<SimulationResult,double> sim_result) const override {
        /* double l1 = 0.0;
        double l2 = 0.0;
        if(sim_result) {
            SimulationResult simResult = *sim_result;
            std::vector<std::pair<SpeciesId,double>>& result = simResult.constrained;
            std::sort(result.begin(), result.end(), comp);
            auto real_conc = this->get_real_concs();
            
            control(result.size() == real_conc.size());
            for (size_t i = 0; i < result.size(); ++i) {
                double a = result[i].second - this->get_real_conc(result[i].first);
                l2 += a*a;

                double b = (simResult.old_values[result[i].first] - result[i].second);
                l1 += b*b;
            }

            const auto& outputs = this->get_outputs();
            
            for (const auto& nc : simResult.not_constrained) {
                if (outputs.find(nc.first) == outputs.end()) {
                    double b = (simResult.old_values[nc.first] - nc.second);
                    l1 += b*b;
                }
            }

            printf("[INFO] L1 = %e\n", p*l1);
            printf("[INFO] L2 = %e\n", (1-p)*l2);

            return S*(p*l1 + (1-p)*l2);

        } else {
            // TODO: scegli l'orizzonte di simulazione
            double error = sim_result.error();
            if(error < 1e-12) {
                error = 1e-12;
            }
            return 1e6*(1/((error)/100.0));
        } */
        TODO("sistemare");

    } 
};

class ErrorSum: public ErrorHandler {

    std::vector<std::pair<ErrorHandler*, double>> handlers;
    double horizon = 100.0;
    double scale = 1e6;

public:
    ~ErrorSum() override = default;

    ErrorSum() = default;

    ErrorSum(double horizon, double scale): horizon(horizon), scale(scale) {}

    void add_handler(ErrorHandler *handler, double scale) {
        this->handlers.push_back(std::pair(handler,scale));
    }

    double error(std::expected<SimulationResult,double> sim_result) const override {
        double result = 0.0;
        if(sim_result) {
            for (const auto& [handler, scale] : this->handlers) {
                result += handler->error(sim_result) * scale;
            }
            return result;
        } else {
            double error = sim_result.error();
            if(error < 1e-12) {
                error = 1e-12;
            }
            return this->scale*(1/((error)/this->horizon));
        }
    }

};

class StabilityError: public ErrorHandler {
public:
    ~StabilityError() override = default;

    StabilityError() = default;

    double error(std::expected<SimulationResult,double> sim_result) const override {
        control(sim_result);
        SimulationResult simResult = *sim_result;

        double result = 0.0;

        auto it_min = simResult.mins.cbegin();
        auto it_max = simResult.maxs.cbegin();
        for (; it_min != simResult.mins.cend() && it_max != simResult.maxs.cend(); ++it_min, ++it_max) {
            double min_val = it_min->second;
            double max_val = it_max->second;

            if (min_val < 0) {
                result += -min_val;
            }
            if (max_val > 1) {
                result += max_val - 1;
            }
        }

        return result;
    }

};

class FittingError: public ErrorHandler {
    std::vector<std::pair<SpeciesId, SpeciesId>> constrain;
    double epsilon = 1e6;
public:
    ~FittingError() override = default;

    FittingError() = default;

    FittingError(double epsilon): epsilon(epsilon) {}

    void add_constrain(SpeciesId id1, SpeciesId id2) {
        constrain.push_back(std::pair(id1,id2));
    }

    double error(std::expected<SimulationResult,double> sim_result) const override {
        control(sim_result);
        SimulationResult simResult = *sim_result;
        double err = 0.0;
        for (const auto &ids : constrain) {
            SpeciesId id1 = ids.first;
            SpeciesId id2 = ids.second;

            auto it1 = simResult.constrained.find(id1);
            auto it2 = simResult.constrained.find(id2);

            if (it1 != simResult.constrained.end() && it2 != simResult.constrained.end()) {
                double diff = it1->second - it2->second;
                double poss_err = epsilon - diff;
                if(poss_err > 0) {
                    err += poss_err;
                }
            } else {
                control(false);
            }
        }
        return err;
    }
};

class TransitorialError: public ErrorHandler {

public:
    ~TransitorialError() override = default;

    TransitorialError() = default;

    double error(std::expected<SimulationResult,double> sim_result) const override {
        control(sim_result);
        SimulationResult simResult = *sim_result;
        double err = 0.0;
        for (const auto& nc : simResult.not_constrained) {
            double b = (simResult.old_values[nc.first] - nc.second);
            err += b*b;
        }
        
        for (const auto& nc : simResult.constrained) {
            double b = (simResult.old_values[nc.first] - nc.second);
            err += b*b;
        }

        return err;
    }

};

#endif // LOSS_FUNCTION_HPP_