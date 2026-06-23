import os
import numpy as np
import matplotlib.pyplot as plt
import traceback
from scipy import io
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

from Neural_Decoding.decoders import LSTMDecoder
from Neural_Decoding.preprocessing_funcs import get_spikes_with_history, bin_spikes

from joblib import Parallel, delayed


def train_and_eval_with_metrics(X_train, y_train, X_test, y_test):
    n_train = min(X_train.shape[0], y_train.shape[0])
    n_test = min(X_test.shape[0], y_test.shape[0])
    X_train = X_train[:n_train]
    y_train = y_train[:n_train]
    X_test = X_test[:n_test]
    y_test = y_test[:n_test]

    model = LSTMDecoder(units=400, num_epochs=5)
    y_train = y_train[:, [0, 1]]
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_test = y_test[:y_pred.shape[0]]
    
    # Filter: compute errors only for samples with speed >= vlim
    valid_idx = y_test[:, 2] >= vlim  # if y_test is (n_samples, k), speed is in column 2
    if not np.any(valid_idx):
        raise ValueError("No evaluation samples remain after speed filtering.")
    y_test_filtered = y_test[:, [0, 1]]
    y_test_filtered = y_test_filtered[valid_idx, :]
    y_pred_filtered = y_pred[valid_idx, :]

    diff = y_test_filtered - y_pred_filtered
    distances = np.linalg.norm(diff, axis=1)
    mse = mean_squared_error(y_test_filtered, y_pred_filtered)
    mae = mean_absolute_error(y_test_filtered, y_pred_filtered)
    r2 = r2_score(y_test_filtered, y_pred_filtered)
    return y_test, y_pred, distances, mse, mae, r2


def cross_test(X1, y1, X2, y2):
    y_real1, y_pred1, dist1, mse1, mae1, r2_1 = train_and_eval_with_metrics(X1, y1, X2, y2)
    y_real2, y_pred2, dist2, mse2, mae2, r2_2 = train_and_eval_with_metrics(X2, y2, X1, y1)
    all_real = np.concatenate([y_real1, y_real2])
    all_predict = np.concatenate([y_pred1, y_pred2])
    all_dists = np.concatenate([dist1, dist2])
    avg_mse = (mse1 + mse2) / 2
    avg_mae = (mae1 + mae2) / 2
    avg_r2  = (r2_1 + r2_2) / 2
    return all_real, all_predict, all_dists, avg_mse, avg_mae, avg_r2


def run_shuffle_iteration(segments, trimmed_rows, neural_data, Set1, Set2,
                          bins_before, bins_after, bins_current, y1, y2):
    segments = segments.copy()
    np.random.shuffle(segments)
    reshuffled = segments.reshape(-1, neural_data.shape[1])
    reshuffled = np.vstack([reshuffled, neural_data[trimmed_rows:]])

    randX1 = get_spikes_with_history(reshuffled[Set1, :], bins_before, bins_after, bins_current)
    randX2 = get_spikes_with_history(reshuffled[Set2, :], bins_before, bins_after, bins_current)

    randX1 = np.nan_to_num(randX1)
    randX2 = np.nan_to_num(randX2)

    return cross_test(randX1, y1, randX2, y2)


def shuffle_with_convergence(segments, trimmed_rows, neural_data,
                              Set1, Set2, bins_before, bins_after, bins_current,
                              y1, y2,
                              max_runs=10, check_interval=5, convergence_thresh=0.01,
                              n_jobs=4):

    dist_means = []
    shuffle_real_all, shuffle_pred_all, shuffle_dists_all, shuffle_mse_list, shuffle_mae_list, shuffle_r2_list = [], [], [], [], [], []

    for start in range(0, max_runs, check_interval):
        # === Run 'check_interval' shuffles in parallel ===
        results = Parallel(n_jobs=n_jobs)(delayed(run_shuffle_iteration)(
            segments.copy(), trimmed_rows, neural_data.copy(),
            Set1, Set2, bins_before, bins_after, bins_current, y1, y2
        ) for _ in range(check_interval))

        for all_real, all_predict, all_dists, avg_mse, avg_mae, avg_r2 in results:
            shuffle_real_all.append(all_real)
            shuffle_pred_all.append(all_predict)
            shuffle_dists_all.append(all_dists)
            shuffle_mse_list.append(avg_mse)
            shuffle_mae_list.append(avg_mae)
            shuffle_r2_list.append(avg_r2)

            current_mean = np.mean(np.concatenate(shuffle_dists_all))
            dist_means.append(current_mean)

        # === Convergence check ===
        if len(dist_means) >= check_interval:
            recent_means = dist_means[-check_interval:]
            change = np.abs(recent_means[-1] - recent_means[0]) / recent_means[0]
            print(f"[Iteration {len(dist_means)}] Mean error: {current_mean:.4f} / Relative change: {change:.4f}")

            if change < convergence_thresh:
                print(f"\n[OK] Converged after {len(dist_means)} iterations (change < {convergence_thresh})")
                break

    return dist_means, shuffle_real_all, shuffle_pred_all, shuffle_dists_all, shuffle_mse_list, shuffle_mae_list, shuffle_r2_list


# === Base folder ===
base_folder = r"path/for/each_animal/data/"
if not os.path.isdir(base_folder):
    raise FileNotFoundError(f"base_folder does not exist: {base_folder}")

# Directory of this script file
current_dir = os.path.dirname(os.path.abspath(__file__))
# Create subfolders for outputs
cov_save_dir = os.path.join(current_dir, "Convergence")
os.makedirs(cov_save_dir, exist_ok=True)  # create if missing

save_dir = os.path.join(current_dir, "Results")
os.makedirs(save_dir, exist_ok=True)      # create if missing


CellNum = 50
vlim = 0  # speed threshold for filtering evaluation samples

for folder_lv1 in os.listdir(base_folder):
    path_lv1 = os.path.join(base_folder, folder_lv1)
    if os.path.isdir(path_lv1):
        for folder_lv2 in os.listdir(path_lv1):
            path_lv2 = os.path.join(path_lv1, folder_lv2)
            if os.path.isdir(path_lv2):

                folder = path_lv2 + os.sep
                parent_name = folder_lv1
                child_name = folder_lv2
                prefix = f"{parent_name}_{child_name}"
                print(f"Processing: {folder}")

                try:
                    data = io.loadmat(folder + 'ST_dF_grid_aut_data.mat')
                    spike_times = np.squeeze(data['lk'])
                    Grid_Cells = np.atleast_1d(np.squeeze(data['Grid_Cells'])).astype(int) - 1

                    if len(spike_times) >= CellNum:

                        for i in range(spike_times.shape[0]):
                            spike_times[i] = np.squeeze(spike_times[i])

                        decoding_vars = io.loadmat(folder + 'decoding_variables.mat')['decoding_variables']
                        vel_times = decoding_vars[:, 0]
                        vels = decoding_vars[:, [1, 2, 3]]  # X, Y, Speed
                        caFR = int(round(float(data['caFr'][0, 0])))
                        dt = 1 / caFR
                        t_start = vel_times[0]
                        t_end = vel_times[-1] + dt * 2
                        Whole_neural_data = bin_spikes(spike_times, dt, t_start, t_end)
                        n_time = min(Whole_neural_data.shape[0], vels.shape[0])
                        Whole_neural_data = Whole_neural_data[:n_time]
                        vels = vels[:n_time]

                        # Decoder history window settings
                        bins_before = caFR   # how many bins before the target time
                        bins_current = 1     # whether to include the concurrent bin
                        bins_after = 0       # how many bins after the target time

                        Set1 = np.arange(0, caFR * 60 * 25)              # first 25 min
                        Set2 = np.arange(caFR * 60 * 25, caFR * 60 * 50) # next 25 min

                        # Sort by normalized spatial information and select the cell set
                        SpatialInfo = io.loadmat(folder + 'Spatial_Info.mat')['Info_sum']

                        all_indices = np.arange(SpatialInfo.shape[0])  # all cell indices (0-based)
                        non_grid_indices = np.setdiff1d(all_indices, Grid_Cells)  # exclude grid cells
                        sorted_non_grid = non_grid_indices[np.argsort(SpatialInfo[non_grid_indices, 2])[::-1]]
                        if sorted_non_grid.size < CellNum:
                            print(f"Skipping: insufficient number of non-grid cells ({sorted_non_grid.size} < {CellNum}).")
                            continue
                        best_cells_indices = sorted_non_grid[:CellNum]
                        neural_data = Whole_neural_data[:, best_cells_indices]
                        if Set2[-1] >= neural_data.shape[0] or Set2[-1] >= vels.shape[0]:
                            print("Skipping: recording is shorter than the requested 50 min analysis window.")
                            continue

                        # Prepare data for shuffling
                        shuffled_neural_data = neural_data.copy()

                        # Define segments for shuffling and shuffle
                        n_segments = 1000
                        total_rows = shuffled_neural_data.shape[0]
                        n_segments = min(n_segments, total_rows)
                        segment_length = total_rows // n_segments
                        if segment_length == 0:
                            print("Skipping: not enough rows for shuffle segmentation.")
                            continue
                        trimmed_rows = segment_length * n_segments
                        trimmed_data = shuffled_neural_data[:trimmed_rows]
                        segments = trimmed_data.reshape(n_segments, segment_length, -1).copy()
                        np.random.shuffle(segments)

                        # Reconstruct
                        shuffled_part = segments.reshape(-1, shuffled_neural_data.shape[1])
                        shuffled_neural_data = np.vstack([shuffled_part, neural_data[trimmed_rows:]])

                        # === Preprocess ===
                        X1 = get_spikes_with_history(neural_data[Set1, :], bins_before, bins_after, bins_current)
                        X2 = get_spikes_with_history(neural_data[Set2, :], bins_before, bins_after, bins_current)
                        randX1 = get_spikes_with_history(shuffled_neural_data[Set1, :], bins_before, bins_after, bins_current)
                        randX2 = get_spikes_with_history(shuffled_neural_data[Set2, :], bins_before, bins_after, bins_current)

                        y1 = vels[Set1, :]
                        y2 = vels[Set2, :]

                        # Handle NaNs
                        X1 = np.nan_to_num(X1)
                        X2 = np.nan_to_num(X2)
                        randX1 = np.nan_to_num(randX1)
                        randX2 = np.nan_to_num(randX2)

                        # === Training & prediction ===
                        if os.path.exists(os.path.join(save_dir, f"{prefix}_lstm_results.npz")):
                            data = np.load(os.path.join(save_dir, f"{prefix}_lstm_results.npz"), allow_pickle=True)
                            
                            normal_real = data['normal_real'].copy()
                            normal_predict = data['normal_predict'].copy()
                            normal_dists = data['normal_dists'].copy()
                            # normal_dists = np.linalg.norm(normal_dists, axis=1)  # patch for Ocean EB2 only
                            normal_mse = data['normal_mse'].copy()
                            normal_mae = data['normal_mae'].copy()
                            normal_r2 = data['normal_r2'].copy()
                            print("normal_dists loaded successfully.")
                        
                        else:
                            print("Start prediction")
                            normal_real, normal_predict, normal_dists, normal_mse, normal_mae, normal_r2 = cross_test(X1, y1, X2, y2)
                            # Save results
                            np.savez(os.path.join(save_dir, f"{prefix}_lstm_results.npz"),
                                normal_real=normal_real,
                                normal_predict=normal_predict,
                                normal_dists=normal_dists,
                                normal_mse=normal_mse,
                                normal_mae=normal_mae,
                                normal_r2=normal_r2)
                        
                        if os.path.exists(os.path.join(save_dir, f"{prefix}_lstm_shuffle_results.npz")):
                            data = np.load(os.path.join(save_dir, f"{prefix}_lstm_shuffle_results.npz"), allow_pickle=True)
                            shuffle_dists_all = data['shuffle_dists_all'].copy()
                            shuffle_all_dists = np.concatenate(shuffle_dists_all)

                            # Filter: compute errors only for samples with speed >= vlim
                            shuffle_real_all = data['shuffle_real_all'].copy()
                            shuffle_all_real = np.concatenate(shuffle_real_all)
                            shuffle_pred_all = data['shuffle_pred_all'].copy()
                            shuffle_all_pred = np.concatenate(shuffle_pred_all)
                            valid_idx = shuffle_all_real[:, 2] >= vlim
                            shuffle_real_filtered = shuffle_all_real[:, [0, 1]]
                            shuffle_real_filtered = shuffle_real_filtered[valid_idx]
                            shuffle_pred_filtered = shuffle_all_pred[valid_idx]

                            diff = shuffle_real_filtered - shuffle_pred_filtered
                            shuffle_all_dists = np.linalg.norm(diff, axis=1)
                            mse = mean_squared_error(shuffle_real_filtered, shuffle_pred_filtered)
                            mae = mean_absolute_error(shuffle_real_filtered, shuffle_pred_filtered)
                            r2 = r2_score(shuffle_real_filtered, shuffle_pred_filtered)
                            print("shuffle_dists_all loaded successfully.")
                            
                        else:
                            # === Shuffle loop with automatic convergence check ===
                            dist_means, shuffle_real_all, shuffle_pred_all, shuffle_dists_all, shuffle_mse_list, shuffle_mae_list, shuffle_r2_list = shuffle_with_convergence(
                                segments, trimmed_rows, neural_data,
                                Set1, Set2, bins_before, bins_after, bins_current,
                                y1, y2,
                                max_runs=25,         # maximum number of shuffles
                                check_interval=5,    # check convergence every N runs
                                convergence_thresh=0.01,  # threshold on relative change of mean error
                                n_jobs=4            # adjust to CPU cores
                            )
                            
                            # Save shuffle results
                            np.savez(os.path.join(save_dir, f"{prefix}_lstm_shuffle_results.npz"),
                                 shuffle_real_all=shuffle_real_all, 
                                 shuffle_pred_all=shuffle_pred_all,
                                 shuffle_dists_all=shuffle_dists_all,
                                 shuffle_mse_list=shuffle_mse_list,
                                 shuffle_mae_list=shuffle_mae_list,
                                 shuffle_r2_list=shuffle_r2_list)
            
                            print(f"Saved: {prefix}_lstm_results.npz / {prefix}_shuffle_convergence.jpg")
                            
                            # Plot & save convergence curve
                            plt.figure(figsize=(8, 6))
                            plt.plot(dist_means, marker='o', color='gray', label='Mean Shuffle Error')
                            plt.axhline(y=dist_means[-1], color='blue', linestyle='--', label='Final Mean')
                            plt.title("Convergence of Shuffle Error")
                            plt.xlabel("Shuffle Iteration")
                            plt.ylabel("Mean Error Distance (cm/s)")
                            plt.grid(True)
                            plt.legend()
                            plt.tight_layout()
                            plt.savefig(os.path.join(cov_save_dir, f"{prefix}_shuffle_convergence.jpg"), dpi=300)
                            plt.close()
            
                            shuffle_all_dists = np.concatenate(shuffle_dists_all)

                        # Plot histograms as probability densities (density=True)
                        plt.figure(figsize=(8, 6))
                        shuffle_all_dists = np.abs(shuffle_all_dists)
                        normal_dists = np.abs(normal_dists)
                        plt.hist(shuffle_all_dists, bins=50, color='lightgray', edgecolor='k', alpha=0.8, label='Shuffled', density=True)
                        plt.hist(normal_dists, bins=50, color='steelblue', edgecolor='k', alpha=0.8, label='Real', density=True)
                        plt.title("LSTM Error Distance Distribution")
                        plt.xlabel("Absolute Error (cm/s)")
                        plt.ylabel("Probability")
                        # plt.xlim(0, 40)  # set X-axis limit if needed
                        plt.legend()
                        plt.grid(True)
                        plt.tight_layout()
                        plt.savefig(os.path.join(save_dir, f"{prefix}_error_distribution.jpg"), dpi=300)
                        plt.close()
           
                    else:
                        print("Skipping: insufficient number of cells.")

                except Exception as e:
                    print(f"Error: {folder} -> {e}")
                    traceback.print_exc()
