package genetic_algorithm

import (
    r "math/rand"
    "log"
    "sync"
    "sort"
    "testing"
)

const (
    genotypeLength      =  5
    populationSizeMin   =  50 
    populationSizeMax   =  100
    genotypeMin         = -10
    genotypeMax         =  10
    kRange              =  5
    mutationProbability =  0.35
    mutationAmplitude   =  10
    tournamentSizeMin   =  2
    winnersSizeMin      =  1
)

func TestGeneratePopulation(t *testing.T) {
    log.Println("Testing GeneratePopulation")
    genotypeTemplate := make([][]int, genotypeLength)
    for populationSize:=populationSizeMin;populationSize<=populationSizeMax;populationSize++ {    
        for i:= range(genotypeTemplate) {
            genotypeTemplate[i] = make([]int, 2)
            genotypeTemplate[i][0] = genotypeMin
            genotypeTemplate[i][1] = genotypeMax
        }
        population := GeneratePopulation(populationSize, genotypeTemplate)
        if len(population) != populationSize {
            t.Errorf("Got a population size of %d when it should have been %d\n", len(population), populationSize)
        }
        for i:= range(population) {
            for j:= range(population[i]) {
                if population[i][j] > genotypeMax {
                    t.Errorf("%d-th gene of %d-th individual was beyond the max of %d with value %d\n", i, j, genotypeMax, population[i][j])
                }
                if population[i][j] < genotypeMin {
                    t.Errorf("%d-th gene of %d-th individual was below the min of %d with value %d\n", i, j, genotypeMin, population[i][j])
                }
            }
        }
    }
}

func TestRouletteRanking(t *testing.T) {
    log.Println("Testing RouletteRanking")
    genotypeTemplate := make([][]int, genotypeLength)
    for i:= range(genotypeTemplate) {
        genotypeTemplate[i] = make([]int, 2)
        genotypeTemplate[i][0] = genotypeMin
        genotypeTemplate[i][1] = genotypeMax
    }
    for populationSize:=populationSizeMin;populationSize<=populationSizeMax;populationSize++ {    
        population := GeneratePopulation(populationSize, genotypeTemplate)
        fitness := make([]float64, populationSize)
        sort.Float64s(fitness)
        for i:= range(fitness) {
            fitness[i] = r.Float64()
        }
        for i:=winnersSizeMin;i<populationSize;i++ {
            winners := RouletteRanking(population, fitness, -5, i)
            if len(winners)!=i {
                t.Errorf("The number of winners is different from that requested, respectively %d and %d (max fitness)\n", len(winners), i)
            }
            for i:=range(winners) {
                if winners[i]<0 || winners[i]>=populationSize {
                    t.Errorf("A winner is beyond the range of the population with index %d (population size was %d)\n", winners[i], populationSize)
                }
            }
            winners = RouletteRanking(population, fitness, fitness[0], i)
            if len(winners)!=i {
                t.Errorf("The number of winners is different from that requested, respectively %d and %d (max fitness)\n", len(winners), i)
            }
            for i:=range(winners) {
                if winners[i]<0 || winners[i]>=populationSize {
                    t.Errorf("A winner is beyond the range of the population with index %d (population size was %d)\n", winners[i], populationSize)
                }
            }
        }
    }
}

func TestLinearRanking(t *testing.T) {
    log.Println("Testing LinearRanking")
    genotypeTemplate := make([][]int, genotypeLength)
    for i:= range(genotypeTemplate) {
        genotypeTemplate[i] = make([]int, 2)
        genotypeTemplate[i][0] = genotypeMin
        genotypeTemplate[i][1] = genotypeMax
    }
    var selectionPressures [kRange]float64
    for i:=range selectionPressures {
            selectionPressures[i] = 1. + r.Float64()
    }
    for populationSize:=populationSizeMin;populationSize<=populationSizeMax;populationSize++ {    
        population := GeneratePopulation(populationSize, genotypeTemplate)
        fitness := make([]float64, populationSize)
        for i:= range(fitness) {
            fitness[i] = r.Float64()
        }
        for i:=winnersSizeMin;i<=populationSize;i++ {
            for k := range selectionPressures {
                winners := LinearRanking(population, fitness, true, selectionPressures[k], i)
                if len(winners)!=i {
                        t.Errorf("The total number of winners is different from requested, %d vs %d\n", len(winners), i)
                }
                for j:=range(winners) {
                        if winners[j]<0 || winners[j]>=populationSize {
                            t.Errorf("A winner is beyond the range of the population with index %d (population size was %d)\n", winners[j], populationSize)
                        }
                }
            }
        }
    }
}

func TestExponentionalRanking(t *testing.T) {
    log.Println("Testing ExponentialRanking")
    genotypeTemplate := make([][]int, genotypeLength)
    for i:= range(genotypeTemplate) {
        genotypeTemplate[i] = make([]int, 2)
        genotypeTemplate[i][0] = genotypeMin
        genotypeTemplate[i][1] = genotypeMax
    }
    var selectionPressures [kRange]float64
        for i:=range selectionPressures {
            selectionPressures[i] = 0.01 + r.Float64()*(9./100.)
    }
    for populationSize:=populationSizeMin;populationSize<=populationSizeMax;populationSize++ {    
        population := GeneratePopulation(populationSize, genotypeTemplate)
        fitness := make([]float64, populationSize)
        for i:= range(fitness) {
            fitness[i] = r.Float64()
        }
        for i:=winnersSizeMin;i<=populationSize;i++ {
            for k := range selectionPressures {
                winners := ExponentialRanking(population, fitness, true, selectionPressures[k], i)
                if len(winners)!=i {
                        t.Errorf("The total number of winners is different from requested, %d vs %d\n", len(winners), i)
                }
                for j:=range(winners) {
                        if winners[j]<0 || winners[j]>=populationSize {
                            t.Errorf("A winner is beyond the range of the population with index %d (population size was %d)\n", winners[i], populationSize)
                        }
                }
            }
        }
    }
}

func TestTournamentRanking(t *testing.T) {
    log.Println("Testing TournamentRanking")
    genotypeTemplate := make([][]int, genotypeLength)
    for i:= range(genotypeTemplate) {
        genotypeTemplate[i] = make([]int, 2)
        genotypeTemplate[i][0] = genotypeMin
        genotypeTemplate[i][1] = genotypeMax
    }
    //TODO Make the for loops into goroutine and let populationSize loop over its own loop
    var mux sync.Mutex
    for populationSize:= populationSizeMin; populationSize <= populationSizeMax; populationSize++ {    
        population := GeneratePopulation(populationSize, genotypeTemplate)
        fitness := make([]float64, populationSize)
        for i:= range(fitness) {
            fitness[i] = r.Float64()
        }
        sort.Float64s(fitness)
        c := make(chan int, populationSize - tournamentSizeMin + 1)
        for i:=tournamentSizeMin;i<=populationSize;i++ {
            go func(populationSize int, i int, mux sync.Mutex, c chan int) {         
                for j:=winnersSizeMin;j<=populationSize;j++ {
                    winners := TournamentRanking(population, fitness, true, i, j)
                    if len(winners)!=j {
                        mux.Lock()
                        t.Errorf("The total number of winners is different from requested, %d vs %d\n", len(winners), j)
                        mux.Unlock()
                    }
                    for k:=range(winners) {
                        if winners[k]<0 || winners[k]>=populationSize {
                            mux.Lock()
                            t.Errorf("A winner is beyond the range of the population with index %d (population size was %d)\n", winners[k], populationSize)
                            mux.Unlock()
                        }
                    }
                    for k:= range(winners) {
                        if winners[k] > populationSize - i {
                            mux.Lock()
                            t.Errorf("The %d-th to last individual won a tournament of size %d, which should never happen (min fitness)\n", populationSize - winners[k], i)
                            mux.Unlock()
                        }
                    }
                    winners = TournamentRanking(population, fitness, false, i, j)
                    if len(winners)!=j {
                        mux.Lock()
                        t.Errorf("The total number of winners is different from requested, %d vs %d\n", len(winners), j)
                        mux.Unlock()
                    }
                    for k:= range(winners) {
                        if winners[k] + 1 < i {
                            mux.Lock()
                            t.Errorf("The %d-th to last individual won a tournament of size %d, which should never happen (max fitness)\n", winners[k], i)
                            mux.Unlock()
                        }
                    }
                }
                c <- 1
            }(populationSize, i, mux, c)
        }
        for i:=tournamentSizeMin;i<=populationSize;i++ {
            <-c
        }
    }
}

func TestTwoPointsCrossover(t *testing.T) {
    log.Println("Testing TwoPointsCrossover")
    genotypeTemplate := make([][]int, genotypeLength)
    for i:= range(genotypeTemplate) {
        genotypeTemplate[i] = make([]int, 2)
        genotypeTemplate[i][0] = genotypeMin
        genotypeTemplate[i][1] = genotypeMax
    }
    population := GeneratePopulation(2, genotypeTemplate)
    child := TwoPointsCrossover(population[0], population[1])
    breakpoint1, breakpoint2 := genotypeLength, genotypeLength 
    for i:= range(child) {
        if(child[i]!=population[0][i]) {
            if(child[i]!=population[1][i]) {
                t.Errorf("Child's %d-th gene, with value %d, doesn't belong to either parent1, with value %d, or parent2, witch value %d\n", i, child[i], population[0][i], population[1][i])
            }
            breakpoint1 = i
            break
        }
    }
    for i:=breakpoint1;i<genotypeLength;i++ {
        if(child[i]!=population[1][i]) {
            if(child[i]!=population[0][i]) {
                    t.Errorf("Child's %d-th gene, with value %d, doesn't belong to either parent1, with value %d, or parent2, witch value %d\n", i, child[i], population[0][i], population[1][i])
            }
            breakpoint2 = i
            break
        }
    }
    for i:=breakpoint2;i<genotypeLength;i++ {
        if(child[i]!=population[0][i]) {
            t.Errorf("Child's %d-th gene, which is in the last section of its genome, should belong to parent1, but they are different values, respectively %d and %d (parent2's gene is %d)\nBreakpoints were %d and %d\n", i, child[i], population[0][i], population[1][i], breakpoint1, breakpoint2)
        }
    }
}

func TestUniformCrossover(t *testing.T) {
    log.Println("Testing UniformCrossover")
    genotypeTemplate := make([][]int, genotypeLength)
    for i:= range(genotypeTemplate) {
        genotypeTemplate[i] = make([]int, 2)
        genotypeTemplate[i][0] = genotypeMin
        genotypeTemplate[i][1] = genotypeMax
    }
    population := GeneratePopulation(2, genotypeTemplate)
    child := UniformCrossover(population[0], population[1])
    for i:= range(child) {
        if(child[i]!=population[0][i] && child[i]!=population[1][i]) {
            t.Errorf("Child's %d-th gene doesn't belong to either parent, values are respectively %d, %d and %d\n", i, child[i], population[0][i], population[1][i])
        }
    }
}

func TestMutate(t *testing.T) {
    log.Println("Testing Mutate")
    genotypeTemplate := make([][]int, genotypeLength)
    for i:= range(genotypeTemplate) {
        genotypeTemplate[i] = make([]int, 2)
        genotypeTemplate[i][0] = genotypeMin
        genotypeTemplate[i][1] = genotypeMax
    }
    population := GeneratePopulation(1, genotypeTemplate)
    individual_copy := make([]int, genotypeLength)
    copied := copy(individual_copy, population[0])
    if copied != genotypeLength {
        t.Errorf("Hard-fault: couldn't copy an individual correctly. It was supposed to be of length %d but was instead %d long\n", genotypeLength, copied)
    }
    Mutate(individual_copy, mutationProbability, mutationAmplitude)
    for i:= range(individual_copy) {
        if(individual_copy[i] > population[0][i] + mutationAmplitude || individual_copy[i] < population[0][i] - mutationAmplitude) {
            t.Errorf("The mutated %d-th gene was beyond the requested mutation range, the original gene had value %d, the amplitude was %d and the mutated gene is %d\n", i, population[0][i], mutationAmplitude, individual_copy[i])
        }
    }
}
