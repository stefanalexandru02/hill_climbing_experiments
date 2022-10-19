#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <tuple>
#include <math.h>

using namespace std;

#define _USE_MATH_DEFINES

#define precision 5
#define max_steps 5000
#define run_forever 0

class Function {
public:
    Function(float (*_function) (vector<float>), tuple <float, float> _range)
    {
        function = _function;
        range = _range;
    }
    float (*function) (vector<float>);
    tuple <float, float> range;
};

class Functions{
public:
    enum FunctionType{
        deJongType,
        schwefelType,
        rastringinType,
        michalewiczType
    };
    
    static Function getFunction(FunctionType type)
    {
        switch (type) {
            case deJongType:
                return Function(deJong, make_tuple(-5.12, 5.12));
            case rastringinType:
                return Function(rastringin, make_tuple(-5.12, 5.12));
            case schwefelType:
                return Function(schwefel, make_tuple(-500,500));
            case michalewiczType:
                return Function(michalewicz, make_tuple(0, M_PI));
            default:
                throw logic_error("Type not supported");
        }
    }

    static float deJong(vector<float> values)
    {
        float s = 0;
        for(const auto &x: values)
        {
            s += pow(x, 2);
        }
        return s;
    }
    
    static float schwefel(vector<float> values)
    {
        float s = 0;
        for(const auto&x: values)
        {
            s -= (-x * sin(sqrt(abs(x))));
        }
        return s;
    }
    
    static float rastringin(vector<float> values)
    {
        float s = 10 * values.size();
        for(const auto &x : values)
        {
            s += (pow(x, 2) - 10 * cos(2 * M_PI * x));
        }
        return s;
    }
    
    static float michalewicz(vector<float> values)
    {
        int m = 10;
        int m2 = 2 * m;
        float s = 0;
        for(int i = 0 ; i < values.size(); i++)
        {
            s -= sin(values[i]) * pow((sin((i + 1) * pow(values[i], 2)/M_PI)), m2);
        }
        return s;
    }
};

class Algorithm {
public:
    Algorithm(Function _func, int _dimensions) {
        func = _func.function;
        func_range = _func.range;
        dimensions = _dimensions;
        interval = get<1>(func_range) - get<0>(func_range);
        dimensionBitSize = ceil(log2(interval * pow(10, precision)));
        binaryVectorLenght = dimensions * dimensionBitSize;
    }
    
    template <typename T>
    void printVector(vector<T> vector) {
        for(const auto& x: vector)
            cout<<x<<",";
        cout<<endl;
    }
    
    void printBest(bool evaluation = true, bool values = true, bool candidate = false)
    {
        if(evaluation)
            cout<<best_evaluation<<endl;
        if(values)
            printVector(best_values);
        if(candidate)
            printVector(best_candidate);
    }
    
    void runAlgo()
    {
        executeStep();
    }
    
protected:
    random_device rand;
    
    float (*func) (vector<float>);
    float interval;
    int dimensions;
    int dimensionBitSize;
    int binaryVectorLenght;
    tuple <float, float> func_range;
    
    int current_step;
    bool candidate_exists = false;
    vector<bool> candidate;
    vector<float> candidate_values;
    float candidate_evaluation;
    
    float best_evaluation = numeric_limits<float>::max();
    vector<float> best_values;
    vector<bool> best_candidate;
           
    virtual void executeStep()
    {
        throw logic_error("Execute step not implemented");
    }
    
    float evaluateCandidateValues(vector<float> values)
    {
        // TODO track evaluations maybe and add them to exploredPoints?
        return func(values);
    }
    
    vector<bool> generateCandidate()
    {
        mt19937 rand_generator(rand());
        vector<bool> candidate;
        candidate.reserve(binaryVectorLenght);
        for(int i = 0 ; i < binaryVectorLenght; i++)
        {
            uniform_int_distribution<int> distribution(1, 10);
            int num = distribution(rand_generator);
            if(num % 2 == 0)
                candidate.push_back(false);
            else
                candidate.push_back(true);
        }
        return candidate;
    }
    
    vector<float> decodeCandidate(vector<bool> candidate)
    {
        vector<float> decodedResult;
        decodedResult.reserve(dimensions);
        int f = 0;
        int s = dimensionBitSize;
        for(int i = 0; i < dimensions; i++)
        {
            vector<bool> range(candidate.begin() + f, candidate.begin() + s);
            decodedResult.push_back(decodeElement(range));
            f = s;
            s = s + dimensionBitSize;
        }
        return decodedResult;
    }
    
    float decodeElement(vector<bool> element)
    {
        int s = 0;
        for(const auto x:element)
            s = (s << 1) | x;
        float v = (float)s / (pow(2, dimensionBitSize) - 1) * interval + get<0>(func_range);
        return v;
    }
};

class BestImprovementHillClimbingAlgorithm : public Algorithm {
public:
    BestImprovementHillClimbingAlgorithm(Function func, int dimensions) :
        Algorithm(func, dimensions) {
    }
    
protected:
    void executeStep() override
    {
        if(current_step >= max_steps && run_forever == 0)
            return;
        current_step += 1;
        if(candidate_exists == false)
        {
            candidate = generateCandidate();
            candidate_values = decodeCandidate(candidate);
            candidate_evaluation = evaluateCandidateValues(candidate_values);
            candidate_exists = true;
        }
        tuple<float, int> bestNeighbour = getBestImprovementNeighbour(candidate);
        if(get<0>(bestNeighbour) < candidate_evaluation)
        {
            int switched_index = get<1>(bestNeighbour);
            candidate[switched_index] = 1 - candidate[switched_index];
            candidate_values = decodeCandidate(candidate);
            candidate_evaluation = evaluateCandidateValues(candidate_values);
        }
        
        if(candidate_evaluation < best_evaluation)
        {
            best_evaluation = candidate_evaluation;
            best_candidate = candidate;
            best_values = candidate_values;
        }
        
        printBest(true, false, false);
        executeStep();
    }
    
private:
    tuple<float, int> getBestImprovementNeighbour(vector<bool> originalCandidate)
    {
        float best = numeric_limits<float>::max();
        int switched_index = -1;
        for(int i = 0; i < originalCandidate.size(); i++)
        {
            originalCandidate[i] = 1 - originalCandidate[i];
            vector<float> n_values = decodeCandidate(originalCandidate);
            float n_evaluation = evaluateCandidateValues(n_values);
            if(n_evaluation < best)
            {
                best = n_evaluation;
                switched_index = i;
            }
            originalCandidate[i] = 1 - originalCandidate[i];
        }
        return make_tuple(best, switched_index);
    }
};

class FirstImprovementHillClimbingAlgorithm : public Algorithm {
public:
    FirstImprovementHillClimbingAlgorithm(Function func, int dimensions) :
        Algorithm(func, dimensions) {
    }
    
protected:
    void executeStep() override
    {
        if(current_step >= max_steps && run_forever == 0)
            return;
        current_step += 1;
        if(candidate_exists == false)
        {
            candidate = generateCandidate();
            candidate_values = decodeCandidate(candidate);
            candidate_evaluation = evaluateCandidateValues(candidate_values);
            candidate_exists = true;
        }
        tuple<float, int> firstImprovedNeighbour = getFirstImprovementNeighbour(candidate);
        if(get<0>(firstImprovedNeighbour) < candidate_evaluation)
        {
            int switched_index = get<1>(firstImprovedNeighbour);
            candidate[switched_index] = 1 - candidate[switched_index];
            candidate_values = decodeCandidate(candidate);
            candidate_evaluation = evaluateCandidateValues(candidate_values);
        }
        
        if(candidate_evaluation < best_evaluation)
        {
            best_evaluation = candidate_evaluation;
            best_candidate = candidate;
            best_values = candidate_values;
        }
        
        printBest(true, false, false);
        executeStep();
    }
    
private:
    tuple<float, int> getFirstImprovementNeighbour(vector<bool> originalCandidate)
    {
        float best = candidate_evaluation;
        int switched_index = -1;
        for(int i = 0; i < originalCandidate.size(); i++)
        {
            originalCandidate[i] = 1 - originalCandidate[i];
            vector<float> n_values = decodeCandidate(originalCandidate);
            float n_evaluation = evaluateCandidateValues(n_values);
            if(n_evaluation < best)
            {
                best = n_evaluation;
                switched_index = i;
                return make_tuple(best, switched_index);
            }
            originalCandidate[i] = 1 - originalCandidate[i];
        }
        return make_tuple(best, switched_index);
    }
};

class WorstImprovementHillClimbingAlgorithm : public Algorithm {
public:
    WorstImprovementHillClimbingAlgorithm(Function func, int dimensions) :
        Algorithm(func, dimensions) {
    }
    
protected:
    void executeStep() override
    {
        if(current_step >= max_steps && run_forever == 0)
            return;
        current_step += 1;
        if(candidate_exists == false)
        {
            candidate = generateCandidate();
            candidate_values = decodeCandidate(candidate);
            candidate_evaluation = evaluateCandidateValues(candidate_values);
            candidate_exists = true;
        }
        tuple<float, int> firstImprovedNeighbour = getWorstImprovementNeighbour(candidate);
        if(get<0>(firstImprovedNeighbour) < candidate_evaluation)
        {
            int switched_index = get<1>(firstImprovedNeighbour);
            candidate[switched_index] = 1 - candidate[switched_index];
            candidate_values = decodeCandidate(candidate);
            candidate_evaluation = evaluateCandidateValues(candidate_values);
        }
        
        if(candidate_evaluation < best_evaluation)
        {
            best_evaluation = candidate_evaluation;
            best_candidate = candidate;
            best_values = candidate_values;
        }
        
        printBest(true, false, false);
        executeStep();
    }
    
private:
    tuple<float, int> getWorstImprovementNeighbour(vector<bool> originalCandidate)
    {
        vector<tuple<float, int>> improvements;
        float best = numeric_limits<float>::max();
        int switched_index = -1;
        for(int i = 0; i < originalCandidate.size(); i++)
        {
            originalCandidate[i] = 1 - originalCandidate[i];
            vector<float> n_values = decodeCandidate(originalCandidate);
            float n_evaluation = evaluateCandidateValues(n_values);
            if(n_evaluation < best)
            {
                best = n_evaluation;
                switched_index = i;
                improvements.push_back(make_tuple(best, switched_index));
            }
            originalCandidate[i] = 1 - originalCandidate[i];
        }
        sort(improvements.begin(), improvements.end());
        // TODO This might cause an issue if improvements is empty. Handle this
        tuple<float, int> worst_improvement = improvements[improvements.size() - 1];
        return make_tuple(get<0>(worst_improvement), get<1>(worst_improvement));
    }
};

class SimulatedAnnealingHillClimbingAlgorithm : public Algorithm {
public:
    SimulatedAnnealingHillClimbingAlgorithm(Function func, int dimensions) :
        Algorithm(func, dimensions) {
    }
    
protected:
    void executeStep() override
    {
        candidate = generateCandidate();
        candidate_values = decodeCandidate(candidate);
        candidate_evaluation = evaluateCandidateValues(candidate_values);
        
        best_evaluation = candidate_evaluation;
        best_candidate = candidate;
        best_values = candidate_values;
        
        //https://www.researchgate.net/publication/225260290_The_Theory_and_Practice_of_Simulated_Annealing
        float T = 100; // maybe have some better value in here
        
        do{
            int MK = getNumberOfRepetitionsForTemperature(current_step, T);
            for(int i = 0; i < MK; i++)
            {
                candidate_values = decodeCandidate(candidate);
                candidate_evaluation = evaluateCandidateValues(candidate_values);
                
                vector<bool> new_candidate = generateRandomNeighbour(candidate, T);
                vector<float> new_candidate_values = decodeCandidate(new_candidate);
                float new_candidate_evaluation = evaluateCandidateValues(new_candidate_values);
                
                // evaluation should be as small as possible (we need to find the global minimum)
                if(new_candidate_evaluation < candidate_evaluation)
                {
                    candidate = new_candidate;
                }
                else if(getZeroOneRandom() < exp(-abs(new_candidate_evaluation - candidate_evaluation))/T)
                {
                    candidate = new_candidate;
                }
                
                cout<<new_candidate_evaluation<<" "<<best_evaluation<<" "<<T<<" "<<(int)(dimensions * 1 / (100 / T))<<endl;
                //printVector(new_candidate_values);
                
                if(candidate_evaluation < best_evaluation)
                {
                    best_evaluation = candidate_evaluation;
                    best_candidate = candidate;
                    best_values = candidate_values;
                }
            }
            T *= 0.99; // maybe have some better value in here
            current_step++;
        }while(current_step < max_steps);
    }
    
private:
    vector<bool> generateRandomNeighbour(vector<bool> originalCandidate, float temperature)
    {
        mt19937 rand_generator(rand());
        uniform_int_distribution<int> distr(0,(int)originalCandidate.size());
        
        int c = (dimensions * 1 / (100 / temperature));
        if(c <= 0) c = 1;
        
        for(int i = 0 ; i < c; i++)
        {
            int random_index = distr(rand_generator);
            originalCandidate[random_index] = 1 - originalCandidate[random_index];
        }

        return originalCandidate;
    }
    
    int getNumberOfRepetitionsForTemperature(int k, int T)
    {
        // Maybe some other function here
        // return 100 * pow(0.99,k);
        return 10000 * T;
    }
    
    float getZeroOneRandom()
    {
        mt19937 rand_generator(rand());
        uniform_int_distribution<int> distr(100);
        return (float)distr(rand_generator) / 100;
    }
};

int main(int argc, const char * argv[]) {
    // Ceva e prost aici... WI ar trebui sa dea tot bine
    WorstImprovementHillClimbingAlgorithm optimizer(Functions::getFunction(Functions::deJongType), 5);
    optimizer.runAlgo();
    
    cout<<"---------------"<<endl;
    optimizer.printBest(true, true, false);
    return 0;
}
