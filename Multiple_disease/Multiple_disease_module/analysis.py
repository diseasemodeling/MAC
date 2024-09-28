import pickle
from utils import output_dir
import numpy as np
import pandas as pd
from pathlib import Path
import re
from scipy.stats import t
import pdb
import matplotlib.pyplot as plt
from scipy.stats import sem

# Constants
FOLDER_NAMES = ["baseline", "elevated"]
FOLDER_ORDER = {name: i for i, name in enumerate(FOLDER_NAMES)}
SHEET_NAMES = ["HETF", "HETM"]
STEPS = {'step': -30, 
         'step_2010': 3, 
         'step_2017': 10, 
         'step_2027': 20, 
         'step_2037': 30, 
         'step_2047': -1}

INTERVENTION_YEAR = 2018
TIME_UNIT = 2
OUTPUT_DIR = Path(output_dir)


# Function to calculate 95% CI
def calculate_95_ci(data):
    data = data.dropna()  # Drop NaN values to avoid issues in calculation
    if len(data) == 0:
        return np.nan, np.nan
    n = len(data)
    mean = np.mean(data)
    se = np.std(data, ddof=1) / np.sqrt(n)  # Manual calculation of standard error
    h = se * t.ppf((1 + 0.95) / 2, n - 1)
    return mean - h, mean + h

def get_data_from_file(filename):
    with open(filename, 'rb') as file:
        return pickle.load(file)

def extract_random_seed(filename):
    match = re.search(r"RS(\d+)_", filename)
    return match.group(1) if match else None
    
def calculate_percentage_change(new_val, old_val):
    try:
        return ((new_val - old_val) / old_val) * 100
    except ZeroDivisionError:
        return np.nan

def format_statistic(median_val, min_val, max_val, is_val=False):
    if is_val:
        if max_val < 0.01:
            return f"{round(median_val, 4)} ({round(min_val, 4)},{round(max_val, 4)})"
        elif max_val < 1:
            return f"{round(median_val, 3)} ({round(min_val, 3)},{round(max_val, 3)})"
        else:
            return f"{round(median_val, 2)} ({round(min_val, 2)},{round(max_val, 2)})"
    else:
        if max_val < 0.01:
            return round(median_val, 4)
        elif max_val < 1:
            return round(median_val, 3)
        else:
            return round(median_val, 2)

def compute_statistics(i, cancer_results, y_data, HPV_infected, total_pop_age, HIV_new_inf, steps):
    def calculate_statistics(data, step):
        return [x[i][step:].sum() for x in data]

    stats_keys = ["cancer_diagnosis", "cancer_diagnosis_new", "cancer_deaths", "cancer_deaths_new", "cancer_cases", "cancer_cases_new"]
    statistics = {}

    if i in [0]:#[0,2]:
        for key in stats_keys:
            hiv_data = calculate_statistics([x[key]["HIV"] for x in cancer_results], steps['step'])
            nonhiv_data = calculate_statistics([x[key]["nonHIV"] for x in cancer_results], steps['step'])
            total_data = calculate_statistics([x[key]["total"] for x in cancer_results], steps['step'])
            if "_new" not in key:
                statistics[f'RR_{key}'] = np.array(hiv_data) / np.array(nonhiv_data)
            statistics[f'{key}_HIV'] = hiv_data
            statistics[f'{key}_nonHIV'] = nonhiv_data
            statistics[f'{key}_total'] = total_data

        # statistics['cancer_cases_age_CD4'] = [x["cancer_cases_age_CD4"][i, steps['step']*TIME_UNIT] for x in cancer_results]  
        
   
    for year, key in zip([steps['step_2010'], steps['step_2017'], steps['step_2027'], steps['step_2037'], steps['step_2047']], 
                         ['RR_HPV_2010', 'RR_HPV_2017', 'RR_HPV_2027', 'RR_HPV_2037', 'RR_HPV_2047']):
        HPV_HIV = [y["HIV"][i][year] for y in y_data]
        HPV_nonHIV = [y["nonHIV"][i][year] for y in y_data]
        statistics[key] = np.array(HPV_HIV) / np.array(HPV_nonHIV)

    for key in ["HIV", "nonHIV", "total"]:
        for year in [steps['step_2010'], steps['step_2017'], steps['step_2027'], steps['step_2037'], steps['step_2047']]:
            data = [y[key][i][year] for y in y_data]
            new_key = f"HPV_prev_{key}{year}"
            statistics[new_key] = data
   
    for year in [steps['step_2010'], steps['step_2017'], steps['step_2027'], steps['step_2037'], steps['step_2047']]:
        data = [pop['HIV'][i][year] / pop['total'][i][year] for pop in total_pop_age]
        new_key = f"HIV_prev_{year}"
        statistics[new_key] = data

    statistics['HIV_new_inf'] = np.asarray([sum(HIV_new_inf[0][steps['step']:])])

    population_keys = ['HIV_pop', 'nonHIV_pop', 'total_pop', 'HIV_HPV_pop', 'nonHIV_HPV_pop', 'total_HPV_pop']
    for key, pop_data in zip(population_keys, [total_pop_age, total_pop_age, total_pop_age, HPV_infected, HPV_infected, HPV_infected]):
        data_type = key.split('_')[0]
        pop = [pop[data_type][i][steps['step_2047']] for pop in pop_data]
        statistics[key] = pop

    
    return statistics

def convert_to_single_value(stats):
    for key, value in stats.items():
        if isinstance(value, list) or isinstance(value, np.ndarray):
            stats[key] = value[0]
    return stats

def generate_statistics(dfs, baseline_stats, folder_name, random_seed, data):
    for risk_ind in range(len(SHEET_NAMES)):
        stats = compute_statistics(risk_ind, data["cancer_results"], data["y_data"], data["HPV_infected"], data["total_pop_age"], data['HIV_new_inf'], STEPS)
        stats = convert_to_single_value(stats)
        df = pd.DataFrame.from_dict(stats, orient='index', columns=[f"{folder_name}_RS{random_seed}"])
        
        sheet = SHEET_NAMES[risk_ind]
        if sheet in dfs:
            dfs[sheet] = pd.concat([dfs[sheet], df], axis=1)
        else:
            dfs[sheet] = df
        
        if folder_name == "baseline" and sheet not in baseline_stats:
            baseline_stats[sheet] = {}
        
        baseline_stats[sheet][random_seed] = stats if folder_name == "baseline" else baseline_stats[sheet][random_seed]

def summarize_statistics(df, is_ci = True):
    summary_dfs = {}
    folder_names = df.columns.str.extract(r'(.*)_RS')[0]

    for folder_name in folder_names.unique():
        folder_df = df[df.columns[folder_names == folder_name]]
        median_df = folder_df.median(axis=1)
        if is_ci:
            ci_df = folder_df.apply(lambda row: pd.Series(calculate_95_ci(row)), axis=1)
            summary_df = pd.DataFrame({f"{folder_name}": folder_df.apply(lambda x: format_statistic(median_df[x.name], ci_df.loc[x.name][0], ci_df.loc[x.name][1], True), axis=1)})
        else:
            min_df = folder_df.min(axis = 1)
            max_df = folder_df.max(axis = 1)
            summary_df = pd.DataFrame({f"{folder_name}": folder_df.apply(lambda x: format_statistic(median_df[x.name], min_df[x.name], max_df[x.name], True), axis=1)})
        summary_dfs[folder_name] = summary_df

    return pd.concat(summary_dfs.values(), axis=1)

def calculate_percentage_change_stats(df, baseline_df):
    percentage_changes = {f"{folder_name}_RS{random_seed}": calculate_percentage_change(df[f"{folder_name}_RS{random_seed}"], baseline_df.loc[random_seed]) 
                          for random_seed in baseline_df.index 
                          for folder_name in FOLDER_NAMES if folder_name != "baseline" and f"{folder_name}_RS{random_seed}" in df.columns}

    df = pd.DataFrame(percentage_changes).replace([np.inf, -np.inf], np.nan)
    return df

def plot_aggregated_boxplots(df, row_name, scenarios, title, filename):
    # Extract the row of interest
    data_row = df.loc[row_name]
    
    # Prepare data for each scenario
    plot_data = []
    labels = []
    
    for scenario in scenarios:
        # Collect all columns for the current scenario
        scenario_data = data_row.filter(like=scenario)
        plot_data.append(scenario_data.values)
        labels.append(scenario)
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.boxplot(plot_data, labels=labels)
    plt.title(title)
    plt.ylabel('Cancer cases')
    plt.xlabel('Scenario')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR/ f"{filename}.png")

def plot_grouped_boxplots(df, rows, row_name, scenarios, title, filename):
    # Prepare data for each scenario
    plot_data = []
    labels = []

    for ind, row in enumerate(rows):
        for scenario in scenarios:
            # Collect all columns for the current scenario
            scenario_data = df.loc[row].filter(like=scenario)
            plot_data.append(scenario_data.values)
            labels.append(f"{row_name[ind]}")

    # Plotting
    plt.figure(figsize=(12, 8))
    plt.boxplot(plot_data, labels=labels)
    plt.title(title)
    plt.ylabel('Percentage Change')
    plt.xlabel('Metric')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR/f"{filename}.png")


def calc_age_cd4_dist(data, step):
    risk_ind = 0
    return [x["cancer_cases_age_CD4"][risk_ind, step] for x in data] 
    
def plot_cumulative_cancer_cases(cancer_cases):
    # Extract unique folder names
    folder_names = sorted(set(col.split('_')[0] for col in cancer_cases.columns), key=lambda x: FOLDER_ORDER[x])
    plt.figure(figsize=(12, 8))
    output_file_path = OUTPUT_DIR / 'total_HIV_VLS.xlsx'
    with pd.ExcelWriter(output_file_path) as writer:
        for folder in folder_names:
            # Select columns related to the current folder
            folder_data = cancer_cases.filter(regex=f'^{folder}_RS')
            # Calculate mean and SEM (Standard Error of the Mean) for confidence intervals
            median_values = folder_data.median(axis=1)
            sem_values = sem(folder_data, axis=1)

            # Calculate confidence intervals (95% CI)
            ci_upper = median_values + 1.96 * sem_values
            ci_lower = median_values - 1.96 * sem_values

            # Plot median values and fill between CI
            plt.plot(median_values, label=f'{folder} Median')
            plt.fill_between(range(len(median_values)), ci_lower, ci_upper, alpha=0.2)
            folder_data.to_excel(writer, sheet_name=f"{folder}")
    years = list(range(INTERVENTION_YEAR, INTERVENTION_YEAR + len(median_values)))
    xticks = range(0, len(median_values), 5)
    xticklabels = years[::5]
    plt.xticks(ticks=xticks, labels=xticklabels)


    plt.xlabel('Year')
    plt.ylabel('Cumulative Cancer Cases Among HIV-positive Population')
    plt.title('Cumulative Cancer Cases with Confidence Intervals')
    plt.legend()
    plt.show()

def calculate_cumulative_num(risk_ind, data, step,  pop_name, name = None):
    if name is not None:
        return np.cumsum(data[name][pop_name][risk_ind][step:])
    else:
        return np.cumsum(data[pop_name][risk_ind][step:])

def generate_cumulative_cancer_cases(cancer_cases, folder_name, random_seed, data):
    risk_ind = 0
    # cumulative_cancer_cases = calculate_cumulative_num(risk_ind, data["cancer_results"][0], STEPS['step'], ', 'HIV', 'cancer_cases_new)
    cumulative_cancer_cases = calculate_cumulative_num(risk_ind, data["total_pop_age"][0], STEPS['step'], 'HIV')
    df = pd.DataFrame(cumulative_cancer_cases, columns=[f"{folder_name}_RS{random_seed}"])
    cancer_cases = pd.concat([cancer_cases, df], axis=1)
    return cancer_cases

# Main processing
is_ci = False
dfs, baseline_stats = {}, {}
cancer_cases = pd.DataFrame()
for folder_name in FOLDER_NAMES:
    random_seed_stats = {}
    folder_path = OUTPUT_DIR / folder_name

    for filename in folder_path.glob("results_*.pkl"):
        data = get_data_from_file(filename)
        random_seed = extract_random_seed(filename.name)
        if random_seed not in random_seed_stats:
            random_seed_stats[random_seed] = {"cancer_results": [], "y_data": [], "HPV_infected": [], "total_pop_age": [], "HIV_new_inf":[],"HIV_new_diag":[]}
        random_seed_stats[random_seed]["cancer_results"].append(data[0])
        random_seed_stats[random_seed]["y_data"].append(data[1])
        random_seed_stats[random_seed]["HPV_infected"].append(data[2])
        random_seed_stats[random_seed]["total_pop_age"].append(data[3])
        random_seed_stats[random_seed]["HIV_new_inf"].append(data[-4])
        # random_seed_stats[random_seed]["HIV_new_diag"].append(data[-3])

    for random_seed, data in random_seed_stats.items():
        generate_statistics(dfs, baseline_stats, folder_name, random_seed, data)
        cancer_cases = generate_cumulative_cancer_cases(cancer_cases, folder_name, random_seed, data)

# plot_cumulative_cancer_cases(cancer_cases)

# Final summary and percentage change calculations
final_summary_df, percentage_change_dfs = {}, {}
for sheet, df in dfs.items():
    final_summary_df[sheet] = summarize_statistics(df, is_ci)

percentage_change_dfs_plot = {}
for sheet, df in dfs.items():
    if sheet in baseline_stats:
        baseline_df = pd.DataFrame.from_dict(baseline_stats[sheet], orient='index')
        df = calculate_percentage_change_stats(df, baseline_df)
        percentage_change_dfs_plot[sheet] =  df
        percentage_change_dfs[sheet] = summarize_statistics(df, is_ci)


output_file_path = OUTPUT_DIR / 'aggregate_results_with_percentage_change_VLS.xlsx'
with pd.ExcelWriter(output_file_path) as writer:
    for sheet, df in final_summary_df.items():
        df.to_excel(writer, sheet_name=f"{sheet}_values")
    for sheet, df in percentage_change_dfs.items():
        df.to_excel(writer, sheet_name=f"{sheet}_percentage_change")

output_file_path = OUTPUT_DIR / 'raw_results_with_percentage_change_VLS.xlsx'
with pd.ExcelWriter(output_file_path) as writer:
    for sheet, df in dfs.items():
        df.to_excel(writer, sheet_name=f"{sheet}_values")
    for sheet, df in percentage_change_dfs_plot.items():
        df.to_excel(writer, sheet_name=f"{sheet}_percentage_change")

# plot_aggregated_boxplots(dfs['HETF'], 'cancer_cases_HIV', FOLDER_NAMES, 
#                          'Cancer Cases Per 100,000 HIV Women (2018-2037)', "relative_boxplot")
# plot_aggregated_boxplots(dfs['HETF'], 'cancer_cases_new_HIV', FOLDER_NAMES, 
#                          'Cancer Cases Among HIV Women (2018-2037)', "absolute_boxplot")
# plot_aggregated_boxplots(percentage_change_dfs_plot['HETF'], 'cancer_cases_HIV', [i for i in FOLDER_NAMES if i != "baseline"], 
#                          'Percentage Change in Cancer Cases Per 100,000 HIV Women (2018-2037) Compared to Baseline', "relative_boxplot_perc_change")
# plot_aggregated_boxplots(percentage_change_dfs_plot['HETF'], 'cancer_cases_new_HIV', [i for i in FOLDER_NAMES if i != "baseline"], 
#                          'Percentage Change in Cancer Cases Among HIV Women (2018-2037) Compared to Baseline', "absolute_boxplot_perc_change")

# plot_grouped_boxplots(percentage_change_dfs_plot['HETF'], ['HIV_new_inf', 'HIV_prev_-1'], 
#                       ['HIV incidence (2018-2037)', 'HIV prevalence (2037)'],
#                         [i for i in FOLDER_NAMES if i != "baseline"], 
#                             'Percentage Change in HIV Incidence and Prevalence Among Women Compared to Baseline', "HIV_perc_change")
    
