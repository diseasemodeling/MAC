from pathlib import Path
import datetime

project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "Data"
today = datetime.date.today().strftime("%b-%d-%Y")
output_dir = project_dir / "Outputs" / today 
output_dir.mkdir(parents=True, exist_ok=True)
    
current_dir = Path(__file__).resolve().parent
config_path = current_dir / 'config.json'
