import pandas as pd
import numpy as np
from .config_reader import ConfigReader
from .utils import data_dir, config_path
import pdb
config = ConfigReader(config_path)

# define key to access to the excel file 
PROP_ACTS_CASUAL_KEY = "prop_protected_acts_casual_name"
PROP_ACTS_MAIN_KEY = "prop_protected_acts_main_name"
PROP_ACTS_ANAL_KEY = "prop_anal_acts_name"
AGE_MIXING_KEY = "age_mixing_name"
RISK_MIXING_KEY = "risk_mixing_name"
CONDOM_EFFICACY_KEY = "condom_efficacy_name"
NUM_ACTS_KEY = "num_acts_name"

INIT_AGE_DIST_KEY = "init_age_dist"
AGE_LB_KEY = "age_lb"
AGE_UB_KEY = "age_ub"

CONTACT_RATE_CHANGE_KEY = "contact_rate_change_name"
PTNR_CHANE_KEY = "ptnr_changes_name"
DEGREE_BIN_MIXING_KEY = "degree_bin_mix_name"
MIN_DEG_KEY = "min_degree"
MAX_DEG_KEY =  "max_degree"
DEG_DIST_BIN_KEY = "degree_dist_bin_all"
DEG_DIST_KEY = "degree_dist"
JUR_KEY = "num_jurisdictions"
def get_num_ptnrs_key(time_unit):
    if time_unit == 2: NUM_PTNRS_KEY = "num_ptnrs_6months" 
    elif time_unit == 1: NUM_PTNRS_KEY = "num_ptnrs_year" 
    elif time_unit == 12: NUM_PTNRS_KEY = "num_ptnrs_month" 
    else: raise("Error in time unit")
    return NUM_PTNRS_KEY
NUM_PTNRS_CASUAL_KEY = "num_ptnrs_month"   
NUM_RISK_GROUP_KEY = "num_risk_group"


class SexualBehavior:
    def __init__(self, time_unit = 2):
        
        NUM_PTNRS_KEY = get_num_ptnrs_key(time_unit)
        self.init_age_dist = config.get(INIT_AGE_DIST_KEY)
        self.num_risk_group =  config.get(NUM_RISK_GROUP_KEY)
        self.num_age_group = len(self.init_age_dist)
        self.age_lb = config.get(AGE_LB_KEY)
        self.age_ub = config.get(AGE_UB_KEY)
        
        self.min_degree = config.get(MIN_DEG_KEY)
        self.max_degree = config.get(MAX_DEG_KEY)
        self.degree_dist_bin_all = self._parse_list2numpy(config.get(DEG_DIST_BIN_KEY))
        self.degree_dist = self._parse_list2numpy(config.get(DEG_DIST_KEY))
        if self.min_degree ==1: self.modify_degree_dist()
        
        self.num_degree_bin = len(self.degree_dist_bin_all)  
        self.non_zero_degree_bin = int(np.ceil(np.log2(self.min_degree) + 1))
        
        excel_path = data_dir/ config.get("behavior_data")
        sexual_data = pd.ExcelFile(excel_path)

        self.prop_protected_acts_casual = self._parse_excel_transpose(sexual_data, config.get(PROP_ACTS_CASUAL_KEY))
        self.prop_protected_acts_main = self._parse_excel_transpose(sexual_data, config.get(PROP_ACTS_MAIN_KEY))
        self.prop_anal_acts = self._parse_excel_transpose(sexual_data, config.get(PROP_ACTS_ANAL_KEY))
        self.age_mixing = self._parse_age_mixing(sexual_data, config.get(AGE_MIXING_KEY))
        self.risk_mixing = self._parse_risk_mixing(sexual_data, config.get(RISK_MIXING_KEY))
        self.degree_bin_mixing = self._parse_degree_mix(sexual_data, config.get(DEGREE_BIN_MIXING_KEY))
        self.condom_efficacy = self._parse_condom_efficacy(sexual_data, config.get(CONDOM_EFFICACY_KEY))
        self.num_acts = self._parse_excel_transpose(sexual_data, config.get(NUM_ACTS_KEY))
        self.num_ptnrs = self._parse_num_ptnrs(sexual_data,  config.get(NUM_PTNRS_KEY))
        self.num_ptnrs_casual = self._parse_num_ptnrs(sexual_data,  config.get(NUM_PTNRS_CASUAL_KEY))
        self.contact_rate_change = self._parse_contact_rate_change(sexual_data, config.get(CONTACT_RATE_CHANGE_KEY))
        
        self.num_jurisdictions = config.get(JUR_KEY)
        self.calc_jurisdiction_mixing()
        self.calc_age_degree_joint_dist()
        self.calc_ptnrship_mixing()
        
    def _parse_list2numpy(self, list):
        return np.array(list)
        
    def _parse_excel_transpose(self, sexual_data, sheet_name):
        return sexual_data.parse(sheet_name=sheet_name, index_col=0).T.values[:self.num_risk_group]
    
    def _parse_excel(self, sexual_data, sheet_name):
        return sexual_data.parse(sheet_name=sheet_name, index_col=0).values[:self.num_risk_group]
    
    def _parse_age_mixing(self, sexual_data, sheet_name):
        return (sexual_data.parse(sheet_name=sheet_name, header=None)
                .values[:int(self.num_risk_group* self.num_age_group)]
                .reshape((self.num_risk_group, self.num_age_group, self.num_age_group)) / 100)
    
    def _parse_degree_mix(self, sexual_data, sheet_name):
        return (sexual_data.parse(sheet_name=sheet_name, header=None)
                .values[:int(self.num_risk_group* self.num_degree_bin)]
                .reshape((self.num_risk_group, self.num_degree_bin, self.num_degree_bin)))
    
    def _parse_risk_mixing(self, sexual_data, sheet_name):
        return sexual_data.parse(sheet_name=sheet_name, header=None).values
        
    
    def _parse_num_ptnrs(self, sexual_data, sheet_name):
        return (sexual_data.parse(sheet_name=sheet_name, header=None)\
                .values[:int(self.num_risk_group* self.num_age_group)]
                .reshape((self.num_risk_group, self.num_age_group, self.num_degree_bin)))
    
    def _parse_contact_rate_change(self, sexual_data, sheet_name):
        return (sexual_data.parse(sheet_name = sheet_name, header = None)
                .values[:int(self.num_risk_group* self.num_age_group)]
                .reshape((self.num_risk_group, self.num_age_group, self.num_degree_bin)))
    
    def _parse_condom_efficacy(self, sexual_data, sheet_name):
        return sexual_data.parse(sheet_name = sheet_name, index_col = None).T.values[:self.num_risk_group]
    
    def modify_degree_dist(self):
        #set 20% of females have 1 lifetime partner; 8% of males have 1 lifetime partner
        deg1 = config.get("deg1")
        for i in range(self.num_risk_group):
            self.degree_dist[i] *= 1 - deg1[i]
            self.degree_dist[i, 1] = deg1[i]
                
    # calculate partnership mixing 
    def calc_ptnrship_mixing(self):
        risk_mix = self.risk_mixing[:self.num_risk_group,:self.num_risk_group,None,None, None, None, None, None]
        age_mix = self.age_mixing[:, None,:,:, None, None, None, None]
        degree_mix = self.degree_bin_mixing[:, None,None, None, :, :, None, None] # degree_ind, degree_num
        jur_mix = self.jurisdiction_mixing[None, None,None, None, None, None :, :]
        mixing = risk_mix * age_mix * degree_mix * jur_mix # risk, risk, age_ind, age_num, degree_ind, degree_num, jur_ind, jur_num
        contact_rate_change = np.copy(self.contact_rate_change) 
        self.num_adjusted_ptnr = contact_rate_change[:, None, :, None, :, None, None, None] * mixing

    ### calculate jurisdiction mixing matrix
    def calc_jurisdiction_mixing(self):
        self.jurisdiction_list = np.arange(self.num_jurisdictions)
        if self.num_jurisdictions > 1:
            self.jurisdiction_mixing = np.ones((self.num_jurisdictions, self.num_jurisdictions))
            const = (1 - self.jurisdiction_mixing_within)/(self.num_jurisdictions - 1)
            if self.num_jurisdictions > 1:
                self.jurisdiction_mixing = self.jurisdiction_mixing * const
            np.fill_diagonal(self.jurisdiction_mixing, self.jurisdiction_mixing_within)
        else:
            self.jurisdiction_mixing_within = 1
            self.jurisdiction_mixing = np.ones((self.num_jurisdictions, self.num_jurisdictions)) * self.jurisdiction_mixing_within
    
    ### calculate joint distribution of age group and degree bin 
    def calc_age_degree_joint_dist(self):
        temp_var = np.tile(self.init_age_dist, (self.num_risk_group, 1))[:,:,None]
        deg_dist_row = self.degree_dist[:self.num_risk_group, :][:, None, :]
        self.degree_dist_by_age_degree = np.multiply(temp_var, deg_dist_row)


# sb = SexualBehavior()