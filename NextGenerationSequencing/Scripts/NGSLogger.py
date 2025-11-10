import time
from datetime import datetime


class NGSLogger:

    def __init__(self, output_dir):
        self.log_file = open(f"{output_dir}NGS_log.txt", "w")
        self.inital_time = time.time()

    def write_log(self, input_text):
        current_time = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")
        print(f"{current_time}: {input_text}")
        self.log_file.write(f"{current_time}: {input_text}\n")

    def total_time(self):
        current_time = time.time()
        current_time_formatted = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")
        print(f"{current_time_formatted}: Run finished after {round(current_time - self.inital_time, 0)} s.")
        self.log_file.write(
            f"{current_time_formatted}: Run finished after {round(current_time - self.inital_time, 0)} s.\n")
        self.log_file.close()
