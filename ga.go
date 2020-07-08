package genetic_algorithm

import (
    r "math/rand"
    "log"
    "sort"
    "time"
)

func init(){
    r.Seed(time.Now().UnixNano())
}

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

//Given a population of n individuals, it returns the **indexes** of the new generation from the initial one
//by organizing n number of tournaments of tournamentSize participants, each tournament won by the
//participant with the highest fitness, as read in its record in the fitness array
func TournamentRanking(population [][]int, fitness []float64, minOrMax bool, tournamentSize int) []int {
    individuals := len(population)
    if tournamentSize < 2 {
        log.Println("The tournament size is lesser than 2, capping it")
        tournamentSize = 2

    }
    if tournamentSize > individuals {
        log.Println("The tournament size is greater than the number of individuals in the population, capping it")
        tournamentSize = individuals
    }
    var fitnessMultiplier float64
    if minOrMax {
        fitnessMultiplier = -1.
    } else {
        fitnessMultiplier = 1.
    }
    winners := make([]int, individuals)
    for i := 0; i < individuals; i++ {
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
    }
    return winners
}

//Returns an individual that has a genome whose first and last part are parent1's, the middle one is parent2's.
//The cutting points are random
func TwoPointsCrossover(parent1 []int, parent2 []int) []int {
    length := len(parent1)
    if length != len(parent2) {
        log.Fatal("Tried to crossover two individuals with genomes of different length")
    }
    child := make([]int, length)
    x := r.Intn(length)
    y := r.Intn(length)
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

//Returns an individual that has a genome where each gene is randomly taken from either parent1 or parent2
func UniformCrossover(parent1 []int, parent2 []int) []int {
    length := len(parent1)
    if length != len(parent2) {
        log.Fatal("Tried to crossover two individuals with genomes of different length")
    }
    child := make([]int, length)
    for i:=0;i<length;i++ {
        maskBit := r.Intn(2)
        child[i] = parent1[i]*maskBit + parent2[i]*(1 - maskBit)
    }
    return child
}

//Mutates the individual **in place**, each gene has a mutationProbability chance of mutating with a mutation amplitude
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
    for i:= range(individual) {
        if r.Float32()<mutationProbability {
            individual[i] += (r.Intn(mutationAmplitude) + r.Intn(mutationAmplitude))/2
        }
    }
}
