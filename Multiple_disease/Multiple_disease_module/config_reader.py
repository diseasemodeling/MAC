# config_reader.py
import json

class ConfigReader:
    def __init__(self, config_file_path):
        with open(config_file_path, 'r') as f:
            self.config_data = json.load(f)

    def get(self, key):
        return self.config_data.get(key)
