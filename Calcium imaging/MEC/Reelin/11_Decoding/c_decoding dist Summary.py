import os
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score


script_dir = os.path.dirname(os.path.abspath(__file__))
result_folder = os.path.join(script_dir, "Results")
output_excel = os.path.join(script_dir, "lstm_noGrid_bestSPi50cell_errorDist_with_stats.xlsx")
speed_threshold = 0


def flatten_saved_array(array, n_cols):
    array = np.asarray(array)
    if array.dtype == object:
        parts = []
        for item in array.ravel():
            item = np.asarray(item)
            if item.size:
                parts.append(item.reshape(-1, n_cols))
        if not parts:
            return np.empty((0, n_cols))
        return np.vstack(parts)
    return array.reshape(-1, n_cols)


def mean_or_nan(values):
    return float(np.mean(values)) if values.size else np.nan


def std_or_nan(values):
    return float(np.std(values)) if values.size else np.nan


def r2_or_nan(real, pred):
    if real.shape[0] < 2 or pred.shape[0] < 2:
        return np.nan
    return float(r2_score(real, pred))


stats_dict = {
    "mean": {},
    "std": {},
    "r2": {},
    "shuffle_mae_mean": {},
    "shuffle_mae_std": {},
    "shuffle_r2_mean": {},
    "real-shuffle_mae": {},
}

if not os.path.isdir(result_folder):
    raise FileNotFoundError(f"Result folder does not exist: {result_folder}")

file_list = sorted(os.listdir(result_folder))
total_files = len(file_list)
processed_count = 0
for filename in file_list:
    processed_count += 1
    print(f"[{processed_count}/{total_files}] Processing: {filename}")

    if not filename.endswith("_lstm_results.npz"):
        continue

    filepath = os.path.join(result_folder, filename)
    samplename = filename.replace("_lstm_results.npz", "")
    shuffle_filename = filename.replace("_lstm_results.npz", "_lstm_shuffle_results.npz")
    shuffle_filepath = os.path.join(result_folder, shuffle_filename)
    if not os.path.exists(shuffle_filepath):
        print(f"Skipping {samplename}: missing {shuffle_filename}")
        continue

    data = np.load(filepath, allow_pickle=True)
    shuffle_data = np.load(shuffle_filepath, allow_pickle=True)

    if "normal_real" not in data or "normal_predict" not in data:
        print(f"Skipping {samplename}: normal_real or normal_predict is missing")
        continue

    normal_real_all = np.asarray(data["normal_real"])
    normal_pred_all = np.asarray(data["normal_predict"])
    n_normal = min(normal_real_all.shape[0], normal_pred_all.shape[0])
    normal_real_all = normal_real_all[:n_normal]
    normal_pred_all = normal_pred_all[:n_normal]

    normal_speed = normal_real_all[:, 2]
    valid_idx = normal_speed >= speed_threshold
    real = normal_real_all[valid_idx][:, [0, 1]]
    pred = normal_pred_all[valid_idx]
    distances_error = np.linalg.norm(np.abs(real - pred), axis=1)

    shuffle_real_all = flatten_saved_array(shuffle_data["shuffle_real_all"], 3)
    shuffle_pred_all = flatten_saved_array(shuffle_data["shuffle_pred_all"], 2)
    n_shuffle = min(shuffle_real_all.shape[0], shuffle_pred_all.shape[0])
    shuffle_real_all = shuffle_real_all[:n_shuffle]
    shuffle_pred_all = shuffle_pred_all[:n_shuffle]

    shuffle_speed = shuffle_real_all[:, 2]
    shuffle_valid = shuffle_speed >= speed_threshold
    shuffle_real = shuffle_real_all[shuffle_valid][:, [0, 1]]
    shuffle_pred = shuffle_pred_all[shuffle_valid]
    distances_shuffle_error = np.linalg.norm(np.abs(shuffle_real - shuffle_pred), axis=1)

    stats_dict["mean"][samplename] = mean_or_nan(distances_error)
    stats_dict["std"][samplename] = std_or_nan(distances_error)
    stats_dict["r2"][samplename] = r2_or_nan(real, pred)
    stats_dict["shuffle_mae_mean"][samplename] = mean_or_nan(distances_shuffle_error)
    stats_dict["shuffle_mae_std"][samplename] = std_or_nan(distances_shuffle_error)
    stats_dict["shuffle_r2_mean"][samplename] = r2_or_nan(shuffle_real, shuffle_pred)
    stats_dict["real-shuffle_mae"][samplename] = (
        stats_dict["mean"][samplename] - stats_dict["shuffle_mae_mean"][samplename]
    )

df_stats = pd.DataFrame(stats_dict)

with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    df_stats.to_excel(writer, sheet_name="normal_dists_stats")

print(f"Saved successfully: {output_excel}")
