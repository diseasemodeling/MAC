import numpy as np
from .config_reader import ConfigReader
from .disease import Disease
from .sexual_behavior import SexualBehavior
from .utils import config_path
from .output import Output
import pdb
import copy

config = ConfigReader(config_path)  

BIRTH_RATE_KEY = "birth_rate"
NUM_NODES_KEY = "num_nodes"
AGE_TRANSITION_RATE_KEY ="age_transition_rate"
NONAGENTS_INDEX_KEY = "nonagents_index"


class Model():
    def __init__(self, time_unit = 2, 
                       mode = "netlogo",
                       increased_risk = 0,
                       initial_year = 2006,
                       dryrun_step = 70,
                       simulation_step = 10,
                       model_screening = False,
                       model_vaccine = False,
                       plot = False,
                       random_seed = 1076):
        self.time_unit = time_unit
        self.dryrun_step = dryrun_step
        self.simulation_step = simulation_step
        self.mode = mode
        self.increased_risk = increased_risk
        self.random_seed = random_seed

        self.plot = plot 
        
        self.model_screening = model_screening
        self.model_vaccine = model_vaccine
        
        self.num_nodes = config.get(NUM_NODES_KEY)
        self.birth_rate = config.get(BIRTH_RATE_KEY)
        self.age_transition_rate = np.asarray(config.get(AGE_TRANSITION_RATE_KEY))
        self.nonagents_state = config.get(NONAGENTS_INDEX_KEY)
        
        self.num_disease_state_base, self.num_disease_state_model = Disease.get_num_states()

        self.setup(initial_year)
        
    
    def go(self, num_births = None):
        # Increment time (6-month period)
        self.time += 1

        # Check if a year has passed (every time units of periods)
        if ((self.time - 1)% self.time_unit == 0):
            # Increment the year counter 
            if self.index >= self.dryrun_step: 
                self.t += 1
                # if self.index > self.dryrun_step:
                #     self.output.model_matrix_time.append(copy.deepcopy(self.model_matrix))
            if num_births is None: 
                num_births = self.num_nodes * self.birth_rate       # calculate number of births
            self.aging_population(num_births)                # aging population 
            self.index += 1 
        # if self.index > self.dryrun_step: 
        #     self.output.model_matrix_time.append(copy.deepcopy(self.model_matrix))
        self.calc_pop_size_age_state()
        self.disease.transmission()

        screening_condition = (self.dryrun_step > 45) and (self.index > 44) and self.model_screening
        vaccination_condition = self.index >= self.dryrun_step and self.model_vaccine
        self.disease.transition(screening_condition, vaccination_condition, year=self.t, random_seed=self.random_seed)
        if self.mode != "netlogo": self.output_variables()


    def setup(self,initial_year):
        self.initialize_constant_variables()
        self.initialize_model_variables(initial_year = initial_year)
        self.disease = Disease(sb = self.sb,
                               model_matrix = self.model_matrix, 
                               model_change = self.model_change,
                               netlogo_change= self.netlogo_change,
                               time_unit= self.time_unit, 
                               nonagents_state = self.nonagents_state, 
                               mode = self.mode,
                               increased_risk= self.increased_risk)
        self.initialize_pop_dist()
        self.output = Output(total_sim = self.total_sim_len, 
                             total_pop=self.total_pop_age_state, 
                             model_matrix=self.model_matrix,
                             num_risk_group=self.num_risk_group)
    
    def run(self, time):
        for _ in range(time):
            self.go()

    def dryrun(self, num_nodes = None):
        self.scale_degree_to_population(num_nodes)
        self.disease.initialize_infection2population()
        self.run(time = self.dryrun_step * self.time_unit)
        
    
    def output_variables(self):
        temp_time = self.time - self.time_unit
        time = temp_time - (self.dryrun_step -1) * self.time_unit - 1
        if self.index > self.dryrun_step:
            for i in self.output.output_risk:
                self.output.new_cancer_cases_age[i][time], \
                self.output.new_cancer_diagnosis_age[i][time], \
                self.output.new_cancer_deaths_age[i][time], \
                self.output.new_natural_deaths_age[i][time] = self.disease.calc_cancer_related_results(risk_ind = i)
            for i in range(self.num_risk_group):    
                self.output.total_infected_age[i][time] = self.disease.calc_total_infected(risk_ind = i)
            if self.index == self.dryrun_step + 1 and temp_time % self.time_unit == 0 and self.plot:
                self.handle_plotting()    
        else: 
            temp_time = self.time             

    def handle_plotting(self):
        genotype_freqs = self.disease.calc_genotype_frequency(risk_ind=0)
        self.output.plot_genotype_frequency(*genotype_freqs)
        step = 0
        risk_ind = 0
        cancer_stats = self.output.calc_cancer_inc_mort(risk_ind, self.time_unit, self.total_pop_age_state, step)
        if self.model_screening:
            self.output.plot_cervical_cancer_postscreening(*cancer_stats)
        else:
            self.output.plot_cervical_cancer_prescreening(*cancer_stats)


    def initialize_constant_variables(self):
        self.sb = SexualBehavior()
        self.num_age_group = self.sb.num_age_group 
        self.num_degree_bin = self.sb.num_degree_bin
        self.num_jurisdictions = self.sb.num_jurisdictions
        self.num_risk_group = self.sb.num_risk_group
        self.age_lb = self.sb.age_lb
        self.age_ub = self.sb.age_ub
        self.degree_dist_bin_all = self.sb.degree_dist_bin_all
        self.degree_dist_by_age_degree = self.sb.degree_dist_by_age_degree
        
        
    def initialize_model_variables(self, initial_year):
        self.index = 0
        self.time = 0
        self.t = initial_year
        self.nonagents_matrix = np.zeros((self.num_risk_group, self.num_age_group, self.num_degree_bin))
        
        self.model_matrix = {n: np.zeros((self.num_disease_state_base[n], 
                                         self.num_age_group, 
                                         self.num_degree_bin, self.num_jurisdictions, 
                                         self.num_disease_state_model[n])) 
                                         for n in range(self.num_risk_group)}
        self.model_change = {n: np.zeros((self.num_disease_state_base[n], 
                                         self.num_age_group, 
                                         self.num_degree_bin, self.num_jurisdictions, 
                                         self.num_disease_state_model[n],
                                         self.num_disease_state_model[n])) 
                                         for n in range(self.num_risk_group)}
        self.netlogo_change = {n: np.zeros((self.num_disease_state_base[n], 
                                         self.num_age_group, 
                                         self.num_degree_bin, self.num_jurisdictions, 
                                         self.num_disease_state_model[n],
                                         self.num_disease_state_model[n])) 
                                         for n in range(self.num_risk_group)}
        self.init_pop_dist ={n: np.zeros_like(self.model_matrix[n]) for n in range(self.num_risk_group)}
        
        self.total_sim_len = self.simulation_step * self.time_unit #(self.simulation_step + 1) * self.time_unit 
        self.total_pop_age_state = {n: np.zeros((self.total_sim_len,
                                                 self.num_disease_state_base[n], 
                                                 self.num_age_group))
                                    for n in range(self.num_risk_group)}
       
    # initialize population distribution 
    def initialize_pop_dist(self):
        expanded_matrix = self.degree_dist_by_age_degree[...,None]/self.num_jurisdictions
        for risk_ind in range(self.num_risk_group): 
            susceptible_ind = self.disease.susceptible_index[risk_ind]
            self.init_pop_dist[risk_ind][self.nonagents_state, ..., susceptible_ind] = expanded_matrix[risk_ind][None, ...]
    
    def check_population(self,script = None):
        for risk_ind in range(self.num_risk_group):
            if np.any(self.model_matrix[risk_ind] <0): 
                print(script, np.where(self.model_matrix[risk_ind]<0))
    
    #### scale population by risk group, age group and degree bin ####
    ### verified against PATH model output
    def scale_degree_to_population(self, num_nodes = None):
        if num_nodes is not None: self.num_nodes =  np.asarray(num_nodes)
        else: self.num_nodes = np.asarray(self.num_nodes)
        num_nodes_expanded = self.num_nodes[:, None,None, None, None, None]
        for risk_ind in range(self.num_risk_group):
            is_all_zero = np.all((self.model_matrix[risk_ind] == 0))
            if is_all_zero:
                self.model_matrix[risk_ind] = np.round(self.init_pop_dist[risk_ind] * num_nodes_expanded[risk_ind]) 
            
        self.update_nonagents_matrix()    
    
    #### aging nonagents population every year #####
    #### verified against PATH model output
    def aging_population(self, num_births):
        for risk_ind in range(self.num_risk_group):
            death_id_model = self.disease.death_state[self.disease.modeling_disease_id][risk_ind]
            ## initialization
            temp_degree_mat = np.delete(self.model_matrix[risk_ind][self.nonagents_state],death_id_model,axis =-1)
            movement_from_age_group = np.zeros_like(temp_degree_mat)
            movement_to_age_group =  np.zeros_like(temp_degree_mat)
            
            ### movement ###
            movement_from_age_group = temp_degree_mat * self.age_transition_rate[:, None, None, None]
            
            birth_index = 0
            age_out_index = -1
            # new birth
            temp_birth = (self.sb.degree_dist[risk_ind] * num_births[risk_ind]).reshape((1, self.sb.num_degree_bin, 1)) / self.num_jurisdictions
            # only goes to susceptible for the new birth
            movement_to_age_group[birth_index,..., self.disease.susceptible_index[risk_ind]] = temp_birth
            movement_to_age_group[birth_index+1:] = np.copy(movement_from_age_group[:age_out_index])

            # record aging out
            # self.age_deaths_out_counter[risk_ind] += movement_from_age_group[age_out_index].sum()
            
            # update model matrix 
            total_movement = movement_to_age_group - movement_from_age_group
            temp_degree_mat += total_movement 
            # in case of negative
            if np.any(temp_degree_mat < 0): 
                # print("negative values in aging movement")
                # self.count_error_negative[self.index] += (temp_degree_mat < 0).sum()
                temp_degree_mat[np.where(temp_degree_mat < 0)] = 0 
            self.model_matrix[risk_ind][self.nonagents_state,...,:death_id_model] =  temp_degree_mat
            
        self.update_nonagents_matrix()
    
    # update suscepitble matrix that is used in the NetLogo
    def update_nonagents_matrix(self):
        # exclude those in the death state
        for risk_ind in range(self.num_risk_group):
            death_id_model = self.disease.death_state[self.disease.modeling_disease_id][risk_ind]
            num_deaths = (self.model_matrix[risk_ind][self.nonagents_state,..., death_id_model]).sum(-1)
            total = self.model_matrix[risk_ind][self.nonagents_state].sum((2,3))
            self.nonagents_matrix[risk_ind] = total - num_deaths
    
    # calculate total alive population by age and HIV state 
    def calc_pop_size_age_state(self):
        temp_time = self.time - self.time_unit
        time = temp_time - (self.dryrun_step -1) * self.time_unit - 1
        if self.index > self.dryrun_step: 
            for risk_ind in range(self.num_risk_group): 
                death_id_model = self.disease.death_state[self.disease.modeling_disease_id][risk_ind]
                temp_model_matrix = np.delete(self.model_matrix[risk_ind],death_id_model,axis =-1)
                self.total_pop_age_state[risk_ind][time] = temp_model_matrix.sum((2,3,4))

