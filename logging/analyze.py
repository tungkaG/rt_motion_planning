import os
import glob
import pandas as pd

# Output file
output_file = "latency_results.txt"

# Open the output file for writing
with open(output_file, "w") as f_out:
    # Find all CSV files in the current directory
    for csv_file in glob.glob("*.csv"):
        f_out.write(f"Results for {csv_file}\n")
        print(f"Processing {csv_file}...")

        # Read CSV into pandas DataFrame
        # Assuming columns: sample_latency, ignore, traj_gen_latency, traj_eval_latency
        df = pd.read_csv(csv_file, header=None,
                         names=["sample_latency", "ignore", "trajectory_generation_latency", "trajectory_evaluation_latency"])

        # Drop the "ignore" column
        df = df.drop(columns=["ignore"])

        # Process each latency column
        for col in df.columns:
            values = df[col]
            min_val = values.min()
            max_val = values.max()
            avg_val = values.mean()
            std_dev = values.std()
            jitter = max_val - min_val  # jitter defined as range

            # Write results
            f_out.write(f"  {col}:\n")
            f_out.write(f"    Min: {min_val:.3f}\n")
            f_out.write(f"    Max: {max_val:.3f}\n")
            f_out.write(f"    Avg: {avg_val:.3f}\n")
            f_out.write(f"    Std Deviation: {std_dev:.3f}\n")
            f_out.write(f"    Jitter: {jitter:.3f}\n")

        f_out.write("\n")

print(f"Results saved to {output_file}")
