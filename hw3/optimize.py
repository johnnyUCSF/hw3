from hw3 import io
from hw3 import smith_waterman
import random
import copy
#########
def byFirst(elem):
    return elem[0]

def calc_total_scores(scoring_matrix,true_pos,true_neg,pen_gap_open,pen_gap_extend):
    ######note all pairs contains files; each file contains pairs
    all_pairs = [true_pos,true_neg]
    #####
    for file in all_pairs:
        for pair in file:
            ###this contains just the scores
            output = smith_waterman.smith_waterman(pair[0],pair[1],scoring_matrix,pen_gap_open,pen_gap_extend)
            total_scores.append([output[2],pair])
    return(total_scores)

def calc_fitness(FPRs,scoring_matrix,true_pos,true_neg,pen_gap_open,pen_gap_extend):
        ######takes sum of TPRs for each; max = 4.0
    ######note all pairs contains files; each file contains pairs
    all_pairs = [true_pos,true_neg]
    #####calculate alignments based on new scoring matrix
    total_scores = []
    for file in all_pairs:
        for pair in file:
            ####calculate scores
            output = smith_waterman.smith_waterman(pair[0],pair[1],scoring_matrix,pen_gap_open,pen_gap_extend)
            total_scores.append([output[2],pair])
    #####
    total_scores.sort(key=byFirst,reverse=True)
    FP_limits = []
    #####adjust values
    for FPR in FPRs:
        FP_limits.append(FPR * len(true_neg))
    ####init 
    TPR_sum = 0.0
    #####
    for FP_limit in FP_limits:
        TP = 0.0
        FP = 0.0
        P = float(len(true_pos))
        N = float(len(true_neg))
        ####sort biggest to smallest
        total_scores.sort(key=byFirst,reverse=True)
        ####scroll through and count TPs and FPs
        if FP_limit == 0:
            for alignment in total_scores:
                if alignment[1] in true_pos:
                    TP += 1
                else:
                    break
        else:
            for alignment in total_scores:
                if FP >= FP_limit:
                    break
                if alignment[1] in true_pos:
                    TP += 1
                else:
                    FP += 1
        ####save the TPR
        TPR_sum += (TP/P)
    return TPR_sum
########

#########
class Particle:
    ############
    self_C = 0.1
    swarm_C = 1.0
    jitter = 0.5 ###this is a percentage jitter, jitters at percent of speed_init
    speed_init = 1.0
    speed_limit = 10.0
    fitness = 0.0
    fitness_prevbest = 0.0
    identities = []
    x_curr = [[]]
    v_curr = [[]]
    x_prevbest = [[]]
    x_neighborbest = [[]]
    ##########initialize x positions by introducing some jitter
    def initialize_x(self):
        ##
        for row in range(len(self.x_curr)):
            for col in range(len(self.x_curr)):
                ####this maintains the symetrical property of matrix
                if col <= row:
                    ####essentially move them by some amount between 0 and speed init
                    sign = random.choice([True, False])
                    if sign == True:
                        self.x_curr[row][col] = copy.deepcopy(self.x_curr[row][col]) + self.speed_init*random.uniform(0, 1)
                    else:
                        self.x_curr[row][col] = copy.deepcopy(self.x_curr[row][col]) - self.speed_init*random.uniform(0, 1)
                    ####make symetrical
                    tmp  = copy.deepcopy(self.x_curr[row][col])
                    self.x_curr[col][row] = copy.deepcopy(tmp)
        return
    ##########give x_coordinates some jitter
    def jitter_x(self):
        for row in range(len(self.x_curr)):
            for col in range(len(self.x_curr)):
                ####this maintains the symetrical property of matrix
                if col <= row:
                    ####essentially move them by some amount between 0 and speed init
                    sign = random.choice([True, False])
                    if sign == True:
                        self.x_curr[row][col] = copy.deepcopy(self.x_curr[row][col]) + self.jitter*(self.speed_init*random.uniform(0, 1))
                    else:
                        self.x_curr[row][col] = copy.deepcopy(self.x_curr[row][col]) - self.jitter*(self.speed_init*random.uniform(0, 1))
                    ####make symetrical
                    tmp  = copy.deepcopy(self.x_curr[row][col])
                    self.x_curr[col][row] = copy.deepcopy(tmp)
        return
    ##########initialize velocity vector
    def initialize_speed(self):
        ##
        self.v_curr = copy.deepcopy(self.x_curr)
        ##
        for row in range(len(self.x_curr)):
            for col in range(len(self.x_curr)):
                ####this maintains the symetrical property of matrix
                if col <= row:
                    sign = random.choice([True, False])
                    if sign == True:
                        self.v_curr[row][col] = self.speed_init*random.uniform(0, 1)
                    else:
                        self.v_curr[row][col] = -(self.speed_init*random.uniform(0, 1))
                    ####make symetrical
                    tmp  = self.v_curr[row][col]
                    self.v_curr[col][row] = tmp
        return
    ############split off first row from scoring matrix
    def format_identities(self):
        self.identities = copy.deepcopy(self.x_curr[0])
        self.x_curr.pop(0)
        return
    ############rebuild scoring matrix
    def blosum(self):
        tmp = copy.deepcopy(self.x_curr)
        tmp.insert(0,copy.deepcopy(self.identities))
        return(tmp)
    #############update position automatically
    def update_position(self):
        ###initialize random numbers
        random2 = random.choice([True, False])
        random3 = random.choice([True, False])
        ##
        term2 = 0
        term3 = 0
        ##
        if random2:
            #####calculate self movement vector
            term2 = subtract_symetrically(self.x_prevbest,self.x_curr,self.self_C)
        ##
        if random3:
            #####calculate best neighbor movement vector
            term3 = subtract_symetrically(self.x_neighborbest,self.x_curr,self.swarm_C)
        ##
        tmp = add_symetrically(term2,term3,1)
        v_new = add_symetrically(self.v_curr,tmp,1)
        ###assign velocity
        self.v_curr = copy.deepcopy(v_new)
        ###limit speed
        self.limit_velocity()
        ###
        x_new = add_symetrically(self.x_curr,v_new,1)
        ####assign
        self.x_curr = copy.deepcopy(x_new)
        ####
        return
    
    ##############limit velocity vector
    def limit_velocity(self):
        for row in range(len(self.v_curr)):
            for col in range(len(self.v_curr)):
                #####
                if abs(self.v_curr[row][col]) > self.speed_limit:
                    self.v_curr[row][col] = copy.deepcopy(self.speed_limit)
        return
      

###############    
    
def add_symetrically(x_curr,v_curr,multiplier):
    if x_curr == 0:
        return(v_curr)
    if v_curr == 0:
        return(x_curr)
    ####add symetrically
    for row in range(len(x_curr)):
        for col in range(len(x_curr)):
            ####this maintains the symetrical property of matrix
            if col <= row:
                ####add
                x_curr[row][col] = (x_curr[row][col]+v_curr[row][col])*multiplier
                ####make symetrical
                tmp  = x_curr[row][col]
                x_curr[col][row] = tmp
    ########
    return(x_curr)

def subtract_symetrically(x_curr,v_curr,multiplier):
    if x_curr == 0:
        return(v_curr)
    if v_curr == 0:
        return(x_curr)
    ####add symetrically
    for row in range(len(x_curr)):
        for col in range(len(x_curr)):
            ####this maintains the symetrical property of matrix
            if col <= row:
                ####add
                x_curr[row][col] = (x_curr[row][col]-v_curr[row][col])*multiplier
                ####make symetrical
                tmp  = x_curr[row][col]
                x_curr[col][row] = tmp
    ########
    return(x_curr)
#########

def particle_swarm(scoring_matrix,FPRs,pen_gap_extend,pen_gap_open,true_pos,true_neg,num_particles,k_max):
    ######initialize
    All_Particles = []
    initial_fitness = calc_fitness(FPRs,scoring_matrix,true_pos,true_neg,pen_gap_open,pen_gap_extend)
    for i in range(num_particles):
        new_particle = Particle()
        ####start from current blosum matrix
        new_particle.x_curr = copy.deepcopy(scoring_matrix)
        ####split off first row of identities to allow for easy matrix operations
        new_particle.format_identities()
        ####save previous best position
        new_particle.x_prevbest = copy.deepcopy(new_particle.x_curr)
        new_particle.fitness_prevbest = copy.deepcopy(initial_fitness)
        ####introduce some jitter
        new_particle.initialize_x()
        ####randomize velocity matrix
        new_particle.initialize_speed()
        ####calc current fitness
        new_particle.fitness = calc_fitness(FPRs,new_particle.blosum(),true_pos,true_neg,pen_gap_open,pen_gap_extend)
        ####save
        All_Particles.append(new_particle)
    ######
    k = 0
    while k < k_max:
        ####new velocity matrix
        new_particle.initialize_speed()
        ####update previous best locations for each particle
        for particle in All_Particles:
            if particle.fitness > particle.fitness_prevbest:
                particle.fitness_prevbest = copy.deepcopy(particle.fitness)
                particle.x_prevbest = copy.deepcopy(particle.x_curr)
        ####find best particle in group
        max_fitness = -1.0
        max_x_best = []
        for particle in All_Particles:
            if particle.fitness > max_fitness:
                max_fitness = copy.deepcopy(particle.fitness)
                max_x_best = copy.deepcopy(particle.x_curr)
        print(k,max_fitness)
        ####memorize neighborbest to each particle
        for particle in All_Particles:
            particle.x_neighborbest = copy.deepcopy(max_x_best)
        ####for every three , change particle speed
        if k % 3 == 0:
            particle.initialize_speed()
        ####update particle positions
        for particle in All_Particles:
            ###particle.jitter_x()
            particle.update_position()
            particle.fitness = calc_fitness(FPRs,particle.blosum(),true_pos,true_neg,pen_gap_open,pen_gap_extend)
            print(particle.fitness)
        ###
        k += 1

    ########return particle with best fitness
    max_fitness = -1.0
    max_x_neighborbest = []
    for particle in All_Particles:
        if particle.fitness > max_fitness:
            print('okay! ',particle.fitness,max_fitness)
            max_fitness = copy.deepcopy(particle.fitness)
            max_x_neighborbest = copy.deepcopy(particle.blosum())
    return(max_fitness,max_x_neighborbest,All_Particles)
########