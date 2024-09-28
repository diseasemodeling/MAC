import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from .config_reader import ConfigReader
from .utils import data_dir, config_path, output_dir


sns.set_style("dark")
import pickle
from pathlib import Path
import pdb
import time

config = ConfigReader(config_path)

def return_modified_actual_data(actual_data):
    actual_data_LB =  [col for col in actual_data.columns if '(LB)' in col]
    actual_data_LB = actual_data[actual_data_LB]
    actual_data_UB =  [col for col in actual_data.columns if '(UB)' in col]
    actual_data_UB = actual_data[actual_data_UB]
    return actual_data_LB, actual_data_UB
    
# preprocessing data
def calc_mid_age(data, add_col = True):
    data.index = round((data['Age_LB'] + data['Age_UB'])/2 + 1)
    data = data.iloc[:,2:]
    if add_col: data.columns = ['Literature ' + i for i in data.columns]
    return data


class Output():
    def __init__(self, total_sim, total_pop, model_matrix, num_risk_group):
        self.total_sim = total_sim
        self.total_pop_age = total_pop
        self.model_matrix =  model_matrix
        self.initialize_parameters(num_risk_group)
        self.read_validation_data()
        self.initialize_output_variables()

    def initialize_parameters(self, num_risk_group):
        # define general parameters
        self.deg_index = ['0', '2^0','2^1','2^2', '2^3', '2^4', '2^5', '2^6', '2^7']
        self.risk_name = ['HETF', 'HETM', 'MSM']
        self.age_name = ['13-17 years', '18-24 years', '25-29 years', '30-34 years', '35-39 yrs', '40-44 years', '45-64 years',"65-100 years"]
        self.num_age_group = len(self.age_name)
        self.degree = [0,1,2,4,8,16,32,64,128]
        self.mid_degree = [np.mean(self.degree[i:i+2]) for i in range(len(self.degree)-1) ]
        self.mid_degree.insert(0, self.degree[0])
        self.mid_degree[2] = self.degree[2]
        self.num_degree_bin = len(self.deg_index)
        self.num_risk_group = num_risk_group
        self.mid_age = [15, 21, 27, 32, 37, 42, 56, 83]
        self.num_disease_state = config.get("num_disease_state")
        self.output_risk = [0] if self.num_risk_group <= 2 else [0,2]


        # define disease parameters 
        self.pop_conversion =1e5
        self.collapsed_disease_name = ['HPV suscepitble', 'HPV infection', 'CIN', 'Cervical cancer', 'Naturally immuned', 'Vaccine']
        self.group_name = ['HIV-susceptible', 'CD4>500',  'CD4>350', 'CD4>200', 'CD4<200', 'HIV-deaths']
        self.type_name = ['HPV 16/18', 'HI5', 'OHR']
        self.num_HIV_status = len(self.group_name)
        self.num_HPV_type = len(self.type_name)
        self.color = ['#6082B6', '#B2BEB5', '#7393B3', '#8A9A5B', '#A9A9A9', '#36454F']
        
    def read_validation_data(self):
        self.cervical_cancer_data = pd.ExcelFile(data_dir/"cervical_cancer_validation.xlsx")
        self.anal_cancer_data = pd.ExcelFile(data_dir/"anal_cancer_calibration.xlsx")
    
    def initialize_output_variables(self):
        self.new_cancer_cases_age = np.zeros((self.num_risk_group, self.total_sim, self.num_HIV_status, self.num_age_group))
        self.new_cancer_diagnosis_age = np.zeros_like(self.new_cancer_cases_age)
        self.new_cancer_deaths_age = np.zeros_like(self.new_cancer_cases_age)
        self.new_natural_deaths_age  = np.zeros_like(self.new_cancer_cases_age)
        self.total_infected_age =  np.zeros_like(self.new_cancer_cases_age)
        self.total_links = {n:np.zeros((self.total_sim,self.num_HIV_status,self.num_disease_state[1][n])) for n in range(self.num_risk_group)}
        self.total_agents = {n:np.zeros((self.total_sim,self.num_HIV_status,self.num_disease_state[1][n])) for n in range(self.num_risk_group)}
        self.model_matrix_time = []
        self.model_change_time = []
        self.new_HIV_inf = []
        self.new_HIV_diag = []
        # self.netlogo_change_time = []
        
    def _get_actual_data(self, data, sheet_name):
        actual_data = data.parse(sheet_name=sheet_name)
        return calc_mid_age(actual_data, False)
    
    def _configure_plot_settings(self, ax):
        ax.set_xlim(13, 70)
        ax.set_xlabel('Age')
        ax.tick_params(axis='both', which='major', labelsize=12)

    def plot_genotype_frequency_general(self, genotype_frequency_general):
        S16_actual = self._get_actual_data(self.cervical_cancer_data, 'S16')
        S16_LB, S16_UB = return_modified_actual_data(S16_actual)
        
        df = pd.DataFrame(genotype_frequency_general.T, index=self.mid_age, columns=['Simulated ' + i for i in self.type_name])
        fig, axes = plt.subplots(self.num_HPV_type, 1, sharex=True, figsize=(8, 6))

        for i in range(self.num_HPV_type):
            df.iloc[:,i].plot(ax = axes[i], c = self.color[-1], legend = True)
            axes[i].set_title(self.type_name[i])
            axes[i].set_ylabel("HPV Prevalence")
            S16_LB.iloc[:,i].plot(linestyle = 'dashed', ax = axes[i], c= self.color[0], legend = True)
            S16_UB.iloc[:,i].plot(linestyle = 'dashed', ax = axes[i], c= self.color[0], legend = True)
            self._configure_plot_settings(axes[i])
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(output_dir/"genotype_frequency_general.png")
        plt.close()
        
    
    def plot_genotype_frequency_normal(self, genotype_frequency_normal):
        S2_actual = self._get_actual_data(self.cervical_cancer_data, 'S2')
        num_figures = 1
        fig, ax = plt.subplots(num_figures)
        ax.plot(self.mid_age , genotype_frequency_normal.sum(0), c= self.color[-1], label = 'Simulated')
        S2_actual.plot(linestyle = 'dotted', ax = ax, c=self.color[0])
        plt.xlabel('Age',fontsize = 15)
        plt.ylabel('HPV prevalence among normal cytoloty',fontsize = 15)
        plt.legend(fontsize = 12)
        plt.xlim(13,70)
        plt.xticks(fontsize =12)
        plt.yticks(fontsize =12)
        plt.savefig(output_dir/ "genotype_frequency_normal.png")
        plt.close()

    def plot_genotype_frequency_precancer(self, genotype_frequency_precancer):    
        S13_actual = self.cervical_cancer_data.parse(sheet_name = 'S13', index_col = 0)
        S14_actual = self.cervical_cancer_data.parse(sheet_name = 'S14', index_col = 0)
        S15_actual = self.cervical_cancer_data.parse(sheet_name = 'S15', index_col = 0)
        num_precancer_type = genotype_frequency_precancer.shape[0]
        fig, ax= plt.subplots(1,num_precancer_type , figsize = (10,6),sharey = True)
        
        summation = genotype_frequency_precancer.sum(-1)
        normalized = summation/ summation.sum(-1)[:,None]*100
        df = pd.DataFrame(normalized)
        df.columns = self.type_name
        df.index =['Simulated'] * self.num_HPV_type
        S13_actual.T.plot(linestyle = '', marker = 'o',ax = ax[0],legend = False, color = self.color[3]) 
        S14_actual.T.plot(linestyle = '', marker = 'o',ax = ax[1],legend = False, c = self.color[3]) 
        S15_actual.T.plot(linestyle = '', marker = 'o',ax = ax[2],legend = False, c = self.color[3]) 

        for i in range(num_precancer_type ):
            df.iloc[i].plot(linestyle = '',marker = 'v',ax = ax[i], legend = False, c = self.color[-1], fontsize = 12)
            ax[i].set_title("CIN {0}".format(i+1), fontsize = 15)
            ax[i].legend(["Literature (LB)", "Literature (UB)",'Simulated'], fontsize = 12)

        ax[0].set_ylabel("Genotype frequency", fontsize = 15)
        plt.yticks(fontsize = 12)
        plt.xticks(fontsize = 12)
        plt.savefig(output_dir/ "genotype_CIN.png")
        plt.close()
        

    def _plot_data(self, ax, data, linestyle, label, color):
        data.plot(ax=ax, linestyle=linestyle, color=color, label=label)

    def _configure_ax(self, ax, ylabel, legend, ylimit=None):
        ax.set_ylabel(ylabel, fontsize=15)
        ax.legend(legend, fontsize=12)
        if ylimit:
            ax.set_ylim(0, ylimit)
        ax.tick_params(labelsize=12)
    
    def _common_plot(self, actual_inc, actual_mort, cancer_inc, cancer_mort, 
                           title, label, save_name, other_actual_inc = None, other_actual_mort = None):
        fig, ax = plt.subplots(2, 1, sharex=True)

        # Incidence Plot
        self._plot_data(ax[0], actual_inc, 'dashed', label[0], self.color[0])
        if other_actual_inc is not None:
            self._plot_data(ax[0], other_actual_inc, 'dashed', label[1], self.color[1])
        self._plot_data(ax[0], pd.DataFrame(cancer_inc, index=self.mid_age), None, label[-1], self.color[-1])
        
        self._configure_ax(ax[0], "Cervical cancer incidence \n per 100,000 women", label, max(actual_inc)+5)

        # Mortality Plot
        self._plot_data(ax[1], actual_mort, 'dashed', label[0], self.color[0])
        if other_actual_mort is not None:
            self._plot_data(ax[1], other_actual_mort, 'dashed', label[1], self.color[1])
        self._plot_data(ax[1], pd.DataFrame(cancer_mort, index=self.mid_age), None, label[-1], self.color[-1])
        
        self._configure_ax(ax[1], "Cervical cancer mortality \n per 100,000 women", label, max(actual_mort)+5)

        plt.xlim(13, 90)
        plt.suptitle(title, fontsize=15, y=1.12)
        plt.xlabel('Age', fontsize=15)
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(save_name)
        plt.close()
    
    def plot_genotype_frequency(self, genotype_frequency_normal,genotype_frequency_CIN, genotype_frequency_general):
        self.plot_genotype_frequency_general(genotype_frequency_general)
        self.plot_genotype_frequency_normal(genotype_frequency_normal)
        self.plot_genotype_frequency_precancer(genotype_frequency_CIN)
    
    def plot_cervical_cancer_prescreening(self, cancer_inc, cancer_mort):
        S11_actual = self._get_actual_data(self.cervical_cancer_data,'S11')
        S12_actual = self._get_actual_data(self.cervical_cancer_data,'S12')

        self._common_plot(S11_actual['Incidence'], S11_actual['Mortality'], cancer_inc, cancer_mort, 
                          'Cervical cancer incidence and mortality \n per 100,000 women, pre-screening period', 
                          ['CTR (1950-1954)', 'CTR (1960-1964)','Simuated'],
                          output_dir / "cervical_cancer_inc_mort_prescreening.png", 
                          S12_actual['Incidence'], S12_actual['Mortality'])

    def plot_cervical_cancer_postscreening(self, cancer_inc, cancer_mort):
        S17_actual = self._get_actual_data(self.cervical_cancer_data,'S17')
        self._common_plot(S17_actual['Incidence'], S17_actual['Mortality'], cancer_inc, cancer_mort, 
                          'Cervical cancer incidence and mortality \n per 100,000 women, post-screening period', 
                          ['SEER 2006','Simuated'],
                          output_dir / "cervical_cancer_inc_mort_postscreening.png")


    def _aggregate_result(self, data, time_unit):
        new_shape = (data.shape[0] // time_unit, time_unit,) + data.shape[1:]
        new_data = data.reshape(new_shape).sum(axis=1)
        return new_data
    
    
    def process_data_type(self, data_type, risk_index, total_pop, cancer_results, time_unit):
        data = self._aggregate_result(
            getattr(self, f"new_{data_type}")[risk_index], 
            time_unit
        )
        group_slices = {'total': slice(None, -1), 'HIV': slice(1, -1), 'nonHIV': 0}

        # Aggregate and assign results for total, HIV, and nonHIV groups
        for group, sl in group_slices.items():
            cancer_results[data_type][group][risk_index] = self.calculate_age_data(data, total_pop, sl, group)

        # Further calculations if data_type contains age
        if "_age" in data_type:
            simple_data_type = data_type.replace("_age", "")
            new_data_type = f"{simple_data_type}_new"
            for group, sl in group_slices.items():
                cancer_results[simple_data_type][group][risk_index], \
                cancer_results[new_data_type][group][risk_index] = self.calculate_group_data(data, total_pop, sl)

    def calculate_age_data(self, data, total_pop, group_slice, group_name):
        # Perform calculations based on group slice
        if group_name == 'nonHIV':
            return data[:, group_slice] / total_pop[:, group_slice] * self.pop_conversion
        else:
            return data[:, group_slice].sum((1)) / total_pop[:, group_slice].sum((1)) * self.pop_conversion

    def calculate_group_data(self, data, total_pop, group_slice):
        # Aggregate age-specific data for total, HIV, and nonHIV groups
        group_data = data[:, group_slice].sum((1,2)) if group_slice != 0 else data[:, group_slice].sum(1)
        age_total_pop = total_pop[:, group_slice].sum((1,2)) if group_slice != 0 else total_pop[:, group_slice].sum(1)
        return  group_data / age_total_pop * self.pop_conversion, group_data
    
    def output_results(self, random_seed,
                             time_unit, 
                             increased_risk, 
                             modeling_screening, 
                             modeling_vaccine, 
                             write = False):
        
        # Consolidate the data types into a single dictionary of dictionaries
        cancer_result_keys = [
            "cancer_diagnosis_age", "cancer_deaths_age", "cancer_cases_age",
            "cancer_diagnosis", "cancer_deaths", "cancer_cases",
            "cancer_diagnosis_new", "cancer_deaths_new", "cancer_cases_new",
        ]
        groups = ["total", "HIV", "nonHIV"]

        # Initialize cancer results with nested dictionaries
        cancer_results = {key: {group: {} for group in groups} for key in cancer_result_keys}

        for risk_index in self.output_risk:
            total_pop = self.total_pop_age[risk_index][::time_unit]
            
            for data_type in cancer_result_keys:
                if hasattr(self, f"new_{data_type}"):
                    # Process each result data type
                    self.process_data_type(data_type, risk_index, total_pop, cancer_results, time_unit)
        
        cancer_results["cancer_diagnosis_age_CD4"] = self.new_cancer_diagnosis_age
        cancer_results["cancer_deaths_age_CD4"] = self.new_cancer_deaths_age
        cancer_results["cancer_cases_age_CD4"] = self.new_cancer_cases_age
        cancer_results["natural_deaths_age_CD4"] = self.new_natural_deaths_age
        
        HPV_prevalence, HPV_infected, total_population = self.calculate_HPV_prevalence(time_unit)
        CD4_population = self.report_CD4count_results(time_unit)
        
        if write:
            self.save_results_to_disk(random_seed,
                                      cancer_results, 
                                      HPV_prevalence, 
                                      HPV_infected, 
                                      total_population, 
                                      CD4_population, 
                                      increased_risk, 
                                      modeling_screening, 
                                      modeling_vaccine, 
                                      self.new_HIV_inf,
                                      self.new_HIV_diag,
                                      self.model_change_time,
                                    #   self.netlogo_change_time,
                                      self.model_matrix_time) 

    def report_CD4count_results(self, time_unit):
        CD4_population = {}
        for risk_ind in range(self.num_risk_group):
            CD4_population[risk_ind] = self.total_pop_age[risk_ind][::time_unit]
        return CD4_population
        
    def calculate_HPV_prevalence(self, time_unit):                
        HPV_prevalence = {"total": {}, "HIV": {}, "nonHIV": {}}
        HPV_infected = {"total": {}, "HIV": {}, "nonHIV": {}}
        total_population = {"total": {}, "HIV": {}, "nonHIV": {}}
        
        for risk_ind in range(self.num_risk_group):
            total_pop = self.total_pop_age[risk_ind][::time_unit]
            total_infected = self.total_infected_age[risk_ind][::time_unit]
            
            total_population['nonHIV'][risk_ind] =  total_pop[:,0].sum((-1))
            total_population['HIV'][risk_ind] = total_pop[:,1:-1].sum((-1,-2))
            total_population['total'][risk_ind] =  total_pop[:,:-1].sum((-1,-2))
        
            HPV_prevalence['nonHIV'][risk_ind] = total_infected[:,0].sum((-1)) / total_population['nonHIV'][risk_ind]
            HPV_prevalence['HIV'][risk_ind] = total_infected[:,1:-1].sum((-1,-2)) / total_population['HIV'][risk_ind]
            HPV_prevalence['total'][risk_ind] = total_infected[:,:-1].sum((-1,-2)) / total_population['total'][risk_ind]

            HPV_infected['nonHIV'][risk_ind] = total_infected[:,0].sum((-1))
            HPV_infected['HIV'][risk_ind] = total_infected[:,1:-1].sum((-1,-2))
            HPV_infected['total'][risk_ind] = total_infected[:,:-1].sum((-1,-2)) 
        
        return HPV_prevalence, HPV_infected, total_population
            
    def calc_cancer_inc_mort(self, risk_ind, time_unit, total_pop_age_state, step):
        total_pop = total_pop_age_state[risk_ind][::time_unit][step]
        new_cancer_diagnosis_age = self.new_cancer_diagnosis_age[risk_ind][:time_unit].sum((0,1))
        new_cancer_diagnosis = new_cancer_diagnosis_age /total_pop.sum(0) * self.pop_conversion
        new_cancer_deaths_age = self.new_cancer_deaths_age[risk_ind][:time_unit].sum((0,1))
        new_cancer_deaths = new_cancer_deaths_age /total_pop.sum(0) * self.pop_conversion
        return new_cancer_diagnosis, new_cancer_deaths
    
    def save_results_to_disk(self,random_seed, *args):
        timestamp = int(time.time() * 1000)
        path = Path(output_dir/ f"results_RS{random_seed}_{timestamp}.pkl")
        while path.exists():
            num += 1
            path = Path(output_dir/ f"results_RS{random_seed}_{timestamp}.pkl")

        with open(path, 'wb') as file:
            pickle.dump(args, file)

        
