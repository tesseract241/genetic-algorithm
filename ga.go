package genetic_algorithm

import (
    "math"
    r "math/rand"
    "log"
    "sort"
    "time"
)

const (
    stored_probabilities    = 5
    tracked_probabilities   = 10
)

//We store tables of the probabilities for linear and exponential ranking, as they are always the same for a given k2/k1
var (
    linear_ranks_probabilities          [][]float64
    linear_ranks_selection_pressure     [][]float64
    exponential_ranks_probabilities     [][]float64
    exponential_ranks_k1                [][]float64
)

func init(){
    r.Seed(time.Now().UnixNano())
    //To speed up starting time we only allocate the first level of the array
    //and will allocate the second as necessary
    linear_ranks_selection_pressure = make([][]float64, 2)
    exponential_ranks_k1            = make([][]float64, 2)
    for i:= range exponential_ranks_k1 {
        linear_ranks_selection_pressure[i] = make([]float64, tracked_probabilities)
        exponential_ranks_k1[i]            = make([]float64, tracked_probabilities)
    }
    linear_ranks_probabilities      = make([][]float64, tracked_probabilities)
    exponential_ranks_probabilities = make([][]float64, tracked_probabilities)
}

//GeneratePopulation returns a population of size individuals, each having a genome of the same length as genotypeTemplate, 
//and minimum and maximum values of each gene taken from genotypeTemplate
func GeneratePopulation(individuals int, genotypeTemplate [][]int) [][]int {
    if individuals <= 0 {
        log.Fatal("The population needs a positive amount of individuals")
    }
    dimension := len(genotypeTemplate)
    if dimension==0 {
        log.Fatal("The genotype can't be of null length")
    }
    if len(genotypeTemplate[0])!=2 {
        log.Fatal("The template for the genotype needs to have a minimum and a maximum")
    }
    population := make([][]int, individuals)
    for i := range population {
        population[i] = make([]int, dimension)
        for j := range population[i] {
            population[i][j] = r.Intn(genotypeTemplate[j][1] - genotypeTemplate[j][0]) + genotypeTemplate[j][0]
        }
    }
    return population
}

//RouletteRanking picks winners based on a roulette wheel whose sectors' width are proportional to the fitness of each individual.
//If using a minimizing fitness, it needs the minimum fitness, otherwise give a positive number for minFitness
//Returns the indexes of the winners
func RouletteRanking(population [][]int, fitness []float64, minFitness float64, winnersSize int) []int {
    individuals := len(population)
    if individuals == 0 {
        log.Fatal("Requested a RouletteRanking on an empty population\n")
    }
    if len(fitness) != individuals {
        log.Fatalf("The number of fitness values is different from the population size, %d vs %d", len(fitness), individuals)
    }
    cumulativeProbabilities := make([]float64, individuals)
    if minFitness >= 0 {
        cumulativeProbabilities[0] = fitness[0]
        for i:=1;i<individuals;i++ {
            cumulativeProbabilities[i] = cumulativeProbabilities[i-1] + fitness[i]
        }
        reciprocal_sum := 1./fitness[individuals-1]
        for i:= range fitness {
            cumulativeProbabilities[i] *= reciprocal_sum
        }
    } else {
        //For problems in which the fitness must be minimized, the modified fitness = 1/(1 + fitness - minFitness(len(population))
        //is used.
        minFitness/=float64(individuals)
        cumulativeProbabilities[0] = 1./(1. + fitness[0] - minFitness)
        for i:=1;i<individuals;i++ {
            cumulativeProbabilities[i] = cumulativeProbabilities[i-1] + 1./(1. + fitness[i] - minFitness)
        }
        reciprocal_sum := 1./cumulativeProbabilities[individuals-1]
        for i:= range cumulativeProbabilities {
            cumulativeProbabilities[i] *= reciprocal_sum
        }
    }
    winners := make([]int, winnersSize)
    for i:= range winners {
        winners[i] = sort.SearchFloat64s(cumulativeProbabilities, r.Float64())
    }
    return winners
}

//CalculateRanks is an Auxiliary function to calculate the index of each ranked individual
//for the ranking selection functions. Exported because it might be useful to the calling program
func CalculateRanks(fitness []float64, minOrMax bool) []int {
    individuals := len(fitness)
    fitnessOrdered := make([]float64, individuals)
    copy(fitnessOrdered, fitness)
    sort.Float64s(fitnessOrdered)
    inverseFitnessMap := make(map[float64]int, individuals)
    for i, fitness := range fitness {
        inverseFitnessMap[fitness] = i
    }
    //ranksLookup stores the index of the ranked individuals in the population array
    ranksLookup := make([]int, individuals)
    if minOrMax {
        for i:= range ranksLookup {
            ranksLookup[i] = inverseFitnessMap[fitnessOrdered[i]]
        }
    } else {
        for i:= range ranksLookup {
            ranksLookup[individuals - 1 - i] = inverseFitnessMap[fitnessOrdered[i]]
        }
    }
    return ranksLookup
}

//calculateLinearRankingProbabilities provides the cumulative probabilites for linear ranking, 
//following the formula P(r) = k1 − r*k2 with k1=selectionPressure/populationSize 
//and k2=selectionPressure/(populationSize*(populationSize-1))
//NOTE that the returned probabilites are not normalized, both for precision and for speed
func calculateLinearRankingProbabilities(selectionPressure float64, populationSize int) []float64 {
    k2 := selectionPressure/float64(populationSize - 1)
    cumulativeProbabilities := make([]float64, populationSize)
    cumulativeProbabilities[0] = selectionPressure 
    for i:=1;i<populationSize-1;i++ {
        cumulativeProbabilities[i] = cumulativeProbabilities[i-1] + selectionPressure - float64(i)*k2
    }
    //The last P(r) is always zero, but k2*(populationSize-1) can end up more than k1 just as much as it can end up less than that,
    //due to precision errors, so we hard-code it.
    if populationSize > 1 {
        cumulativeProbabilities[populationSize-1] = cumulativeProbabilities[populationSize-2]
    }
    return cumulativeProbabilities
}


//linearRankingProbabilitiesGenerator does the heavy-lifting of checking if we've already got the requested probabilities stored, 
//stores them if we've got empty room for them left, and moves things around to keep the tables
//organized by how much they've been used
func linearRankingProbabilitiesGenerator(selectionPressure float64, populationSize int) []float64 {
    index := -1
    usage :=  1.
    var temp_probabilities []float64
    for i:= range linear_ranks_selection_pressure[0] {
        //NOTE how the condition is different here compared to exponentialRankingProbabilitiesGenerator, because
        //linearRanking changes its value not only based on selectionPressure, but also populationSize, so we need to check both
        if linear_ranks_selection_pressure[0][i]==selectionPressure && len(linear_ranks_probabilities[i])==populationSize{
            index = i
            linear_ranks_selection_pressure[1][i]+=1
            usage = linear_ranks_selection_pressure[1][i]
            break
        }
    }
    hasToSwap := false
    if index>-1 {
        //The selectionPressure is already at least tracked
        if index < stored_probabilities {
            //The selectionPressure is already stored, we check if it's long enough and if we need to sort the tables 
            sorted := true
            if index > 0 && linear_ranks_selection_pressure[1][index-1]<linear_ranks_selection_pressure[1][index] {
                sorted = false
            }
            if sorted {
                //No need to sort them, so we just return the correct table
                return linear_ranks_probabilities[index]
            } else {
                //A sort is needed, we store the wanted table in temp_probabilities
                temp_probabilities = linear_ranks_probabilities[index]
                hasToSwap = true
            }
        } else {
            //The selectionPressure is **not** stored, we need to generate it, we check if its usage has surpassed that
            //of the lowest stored table, and check it for sorting if so
            if usage > linear_ranks_selection_pressure[1][stored_probabilities] {
                temp_probabilities = calculateLinearRankingProbabilities(selectionPressure, populationSize)
                hasToSwap = true
            } else {
                return calculateLinearRankingProbabilities(selectionPressure, populationSize)
            }
        }
    } else {
        //The selectionPressure is not tracked, we check if there are empty tables to put it into, 
        //otherwise we simply generate it and return it
        index = tracked_probabilities
        for i:= range linear_ranks_selection_pressure[1] {
            //We find the last empty table, 
            if linear_ranks_selection_pressure[1][i]==0. {
                index = i
                break
            }
        }
        if index < tracked_probabilities {
            //There's room to at least track it, we add the tracking information
            linear_ranks_selection_pressure[0][index] = selectionPressure
            linear_ranks_selection_pressure[1][index] = 1.
            if index < stored_probabilities {
                //There's also room to store it, we allocate memory and store the table after generating it
                linear_ranks_probabilities[index] = make([]float64, populationSize)
                linear_ranks_probabilities[index] = calculateLinearRankingProbabilities(selectionPressure, populationSize)
                return linear_ranks_probabilities[index]
            }
            //No room to store it, we just return the table without storing it
            return calculateLinearRankingProbabilities(selectionPressure, populationSize)
        } else {
            //No room to track it, we just return the table 
            return calculateLinearRankingProbabilities(selectionPressure, populationSize)
        }
    }
    //If we've identified that the tables need sorting, we keep swapping the table we will return with the one before it
    //until sort.Float64sAreSorted(linear_ranks_selection_pressure[1]) will return true
    for ;hasToSwap;index-- {
        temp_pressure := linear_ranks_selection_pressure[0][index]
        linear_ranks_selection_pressure[0][index] = linear_ranks_selection_pressure[0][index-1]
        linear_ranks_selection_pressure[0][index-1] = temp_pressure
        linear_ranks_selection_pressure[1][index] = linear_ranks_selection_pressure[1][index-1]  
        linear_ranks_selection_pressure[1][index-1] = usage
        linear_ranks_probabilities[index] = linear_ranks_probabilities[index-1]
        linear_ranks_probabilities[index-1] = temp_probabilities
        hasToSwap = index > 1 && linear_ranks_selection_pressure[1][index-2] < linear_ranks_selection_pressure[1][index-1]
    }
    return temp_probabilities
}


//LinearRanking picks winners based on a roulette wheel whose sectors' width is based on the ranks of the individuals, 
//with selectionPressure acting as a parameter to determine how much rank weighs.
//Returns the **indexes** of the winners
func LinearRanking(population [][]int, fitness []float64, minOrMax bool, selectionPressure float64, winnersSize int) []int {
    individuals := len(population)
    if individuals == 0 {
        log.Fatal("Requested a LinearRanking on an empty population\n")
    }
    if len(fitness) != individuals {
        log.Fatalf("The number of fitness values is different from the population size, %d vs %d", len(fitness), individuals)
    }
    if selectionPressure < 1 || selectionPressure > 2 {
        log.Fatalf("The selection pressure must be within 1 and 2, got %f\n", selectionPressure)
    }
    ranksLookup := CalculateRanks(fitness, minOrMax)
    cumulativeProbabilities := linearRankingProbabilitiesGenerator(selectionPressure, individuals)
    winners := make([]int, winnersSize)
    for i:= range winners {
        winners[i] = 
        ranksLookup[sort.SearchFloat64s(cumulativeProbabilities, r.Float64()*cumulativeProbabilities[individuals - 1])]
    }
    return winners
}

//calculateExponentialRankingProbabilities provides the cumulative probabilites for exponential ranking, 
//following the formula P(r) = k1*(1-k1)^r with k1 given
//NOTE that the returned probabilites are not normalized, both for precision and for speed
func calculateExponentialRankingProbabilities(k1 float64, populationSize int, startingValue float64, startingIndex int) []float64 {
    if startingIndex < 0 {
        log.Fatalf("Got a startingIndex of %d, should be non-negative\n", startingIndex)
    }
    if startingValue < 0 {
        log.Fatalf("Got a startingValue of %f, should be non-negative\n", startingValue)
    }
    cumulativeProbabilities := make([]float64, populationSize - startingIndex)
    if startingIndex != 0 {
        cumulativeProbabilities[0] = startingValue + k1*math.Pow(1.-k1, float64(startingIndex+1))
    } else {
        cumulativeProbabilities[0] = k1
    }
    for i:=1;i<populationSize - startingIndex;i++ {
        cumulativeProbabilities[i] = cumulativeProbabilities[i-1] + k1*math.Pow(1.-k1, float64(i+startingIndex))
    }
    return cumulativeProbabilities
}

//exponentialRankingProbabilitiesGenerator does the heavy-lifting of checking if we've already got the requested probabilities stored, 
//stores them if we've got empty room for them left, and moves things around to keep the tables
//organized by how much they've been used
func exponentialRankingProbabilitiesGenerator(k1 float64, populationSize int) []float64 {
    index := -1
    usage :=  1.
    var temp_probabilities []float64
    for i:= range exponential_ranks_k1[0] {
        if exponential_ranks_k1[0][i]==k1 {
            index = i
            exponential_ranks_k1[1][i]+=1
            usage = exponential_ranks_k1[1][i]
            break
        }
    }
    hasToSwap := false
    if index>-1 {
        //The k1 is already at least tracked
        if index < stored_probabilities {
            //The k1 is already stored, we check if we need to sort the tables 
            if len(exponential_ranks_probabilities[index])<populationSize {
                //The stored table is not long enough to cover the populationSize, we extend it by only calculating the new elements
                exponential_ranks_probabilities[index] = append(
                    exponential_ranks_probabilities[index],
                    calculateExponentialRankingProbabilities(
                        k1, populationSize,
                        exponential_ranks_probabilities[index][len(exponential_ranks_probabilities[index])-1],
                        len(exponential_ranks_probabilities[index])-1)...
                    )
            }
            sorted := true
            if index > 0 && exponential_ranks_k1[1][index-1]<exponential_ranks_k1[1][index] {
                sorted = false
            }
            if sorted {
                //No need to sort them, so we just return the correct table
                return exponential_ranks_probabilities[index]
            } else {
                //A sort is needed, we store the wanted table in temp_probabilities
                temp_probabilities = exponential_ranks_probabilities[index]
                hasToSwap = true
            }
        } else {
            //The k1 is **not** stored, we need to generate it, we check if its usage has surpassed that
            //of the lowest table, and check it for sorting if so
            if usage > exponential_ranks_k1[1][stored_probabilities] {
                temp_probabilities = calculateExponentialRankingProbabilities(k1, populationSize, 0, 0)
                hasToSwap = true
            } else {
                return calculateExponentialRankingProbabilities(k1, populationSize, 0, 0)
            }
        }
    } else {
        //The k1 is not tracked, we check if there are empty tables to put it into, 
        //otherwise we simply generate it and return it
        index = tracked_probabilities
        for i:= range exponential_ranks_k1[1] {
            //We find the last empty table, 
            if exponential_ranks_k1[1][i]==0. {
                index = i
                break
            }
        }
        if index < tracked_probabilities {
            //There's room to at least track it, we add the tracking information
            exponential_ranks_k1[0][index] = k1 
            exponential_ranks_k1[1][index] = 1.
            if index < stored_probabilities {
                //There's also room to store it, we allocate memory and store the table after generating it
                exponential_ranks_probabilities[index] = make([]float64, populationSize)
                exponential_ranks_probabilities[index] = calculateExponentialRankingProbabilities(k1, populationSize, 0, 0)
                return exponential_ranks_probabilities[index]
            }
            //No room to store it, we just return the table without storing it
            return calculateExponentialRankingProbabilities(k1, populationSize, 0, 0)
        } else {
            //No room to track it, we just return the table 
            return calculateExponentialRankingProbabilities(k1, populationSize, 0, 0)
        }
    }
    //If we've identified that the tables need sorting, we keep swapping the table we will return with the one before it
    //until sort.Float64sAreSorted(exponential_ranks_k1[1]) will return true
    for ;hasToSwap;index-- {
        temp_pressure := exponential_ranks_k1[0][index]
        exponential_ranks_k1[0][index] = exponential_ranks_k1[0][index-1]
        exponential_ranks_k1[0][index-1] = temp_pressure
        exponential_ranks_k1[1][index] = exponential_ranks_k1[1][index-1]  
        exponential_ranks_k1[1][index-1] = usage
        exponential_ranks_probabilities[index] = exponential_ranks_probabilities[index-1]
        exponential_ranks_probabilities[index-1] = temp_probabilities
        hasToSwap = index > 1 && exponential_ranks_k1[1][index-2] < exponential_ranks_k1[1][index-1]
    }
    return temp_probabilities
}

//ExponentialRanking picks winners based on a roulette wheel whose sectors' width is based on the ranks of the individuals, with k1
//acting as a parameter to determine how much rank weighs. The weight decreases exponentially.
//Returns the **indexes** of the winners
func ExponentialRanking(population [][]int, fitness []float64, minOrMax bool, k1 float64, winnersSize int) []int {
    individuals := len(population)
    if individuals == 0 {
        log.Fatal("Requested a RouletteRanking on an empty population\n")
    }
    if len(fitness) != individuals {
        log.Fatalf("The number of fitness values is different from the population size, %d vs %d", len(fitness), individuals)
    }
    if k1 < 0.01 || k1 > 0.1 {
        log.Fatalf("The k1 must be within 0.01 and 0.1, got %f\n", k1)
    }
    ranksLookup := CalculateRanks(fitness, minOrMax)
    cumulativeProbabilities := exponentialRankingProbabilitiesGenerator(k1, individuals)
    winners := make([]int, winnersSize)
    for i:= range winners {
        winners[i] = 
            ranksLookup[sort.SearchFloat64s(cumulativeProbabilities, r.Float64()*cumulativeProbabilities[individuals - 1])]
    }
    return winners
}

//TournamentRanking returns the **indexes** of the new generation from the initial one of size n
//by organizing n number of tournaments of tournamentSize participants, each tournament won by the
//participant with the highest fitness, as read in its record in the fitness array
//TODO Add sigma-scaling, once you know how it works
func TournamentRanking(population [][]int, fitness []float64, minOrMax bool, tournamentSize int, winnersSize int) []int {
    individuals := len(population)
    if tournamentSize < 2 {
        log.Println("The tournament size is lesser than 2, capping it")
        tournamentSize = 2

    }
    if tournamentSize > individuals {
        log.Println("The tournament size is greater than the number of individuals in the population, capping it")
        tournamentSize = individuals
    }
    ranksLookup := CalculateRanks(fitness, minOrMax)
    winners := make([]int, winnersSize)
    for i:= range winners {
        contenders := r.Perm(individuals)
        winners[i] = contenders[0]
        for j:=1; j<tournamentSize; j++ {
            //TODO Find for what values of individuals/tournamentSize it's worth doing something like
            /*
                if winners[i] < tournamentSize - j {
                    break
                }
            */
            if contenders[j] < winners[i] {
                winners[i] = contenders[j]
            }
        }
    }
    for i:= range winners {
        winners[i] = ranksLookup[winners[i]]
    }
    /*
    var fitnessMultiplier float64
    if minOrMax {
        fitnessMultiplier = -1.
    } else {
        fitnessMultiplier = 1.
    }
    winners := make([]int, winnersSize)
    for i := 0; i < winnersSize; i++ {
        contenders := r.Perm(individuals)
        winner := contenders[0]
        winnerFitness := fitness[winner]
        for j:=1; j < tournamentSize; j++ {
            if fitness[contenders[j]] * fitnessMultiplier > winnerFitness * fitnessMultiplier{
                winner = contenders[j]
                winnerFitness = fitness[contenders[j]]
            }
        }
        winners[i] = winner
    }*/
    return winners
}

//TwoPointsCrossover returns an individual that has a genome whose first and last part are parent1's, the middle one is parent2's.
//The cutting points are random
func TwoPointsCrossover(parent1 []int, parent2 []int) []int {
    length := len(parent1)
    if length != len(parent2) {
        log.Fatal("Tried to crossover two individuals with genomes of different length")
    }
    child := make([]int, length)
    x := r.Intn(length + 1)
    y := r.Intn(length + 1)
    dummy := []int{x, y}
    sort.Ints(dummy)
    min, max := dummy[0], dummy[1]
    for i:=0; i<min; i++ {
        child[i] = parent1[i]
    }
    for i:=min; i<max; i++ {
        child[i] = parent2[i]
    }
    for i:=max; i<length; i++ {
        child[i] = parent1[i]
    }
    return child
}

//UniformCrossover returns an individual that has a genome where each gene is randomly taken from either parent1 or parent2
func UniformCrossover(parent1 []int, parent2 []int) []int {
    length := len(parent1)
    if length != len(parent2) {
        log.Fatal("Tried to crossover two individuals with genomes of different length")
    }
    child := make([]int, length)
    for i:=0; i<length; i++ {
        maskBit := r.Intn(2)
        child[i] = parent1[i]*maskBit + parent2[i]*(1 - maskBit)
    }
    return child
}

//Mutate Mutates the individual **in place**, each gene has a mutationProbability chance of mutating with a mutation amplitude
//triangularly distributed with amplitude of mutationAmplitude
func Mutate(individual []int, mutationProbability float32, mutationAmplitude int) {
    if mutationProbability<0 {
        log.Fatal("The mutation probability should be a positive number")
    }
    if mutationAmplitude == 0 {
        log.Fatal("The mutation amplitude should not be 0")
    }
    if mutationAmplitude < 0 {
        log.Println("The mutation amplitude should be positive, took its absolute value")
    }
    for i:= range individual {
        if r.Float32()<mutationProbability {
            individual[i] += (r.Intn(mutationAmplitude) + r.Intn(mutationAmplitude))/2
        }
    }
}
