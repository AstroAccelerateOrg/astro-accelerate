#ifndef ASTRO_ACCELERATE_AA_PULSCAN_RUNNER_HPP
#define ASTRO_ACCELERATE_AA_PULSCAN_RUNNER_HPP

#include "aa_params.hpp"

#if AA_ENABLE_PULSCAN

#include <cuda_runtime.h>
#include <cufft.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <cerrno>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "aa_device_pulscan.hpp"
#include "aa_gpu_timer.hpp"
#include "aa_log.hpp"
#include "aa_periodicity_strategy.hpp"
#include "aa_timelog.hpp"

#ifndef PIPELINE_ERROR_NO_ERROR
#define PIPELINE_ERROR_NO_ERROR 0
#endif

#ifndef PIPELINE_ERROR_GENERAL_GPU_ERROR
#define PIPELINE_ERROR_GENERAL_GPU_ERROR 5
#endif

#ifndef PIPELINE_ERROR_COPY_TO_HOST
#define PIPELINE_ERROR_COPY_TO_HOST 11
#endif

#ifndef PIPELINE_ERROR_COPY_TO_DEVICE
#define PIPELINE_ERROR_COPY_TO_DEVICE 12
#endif

namespace astroaccelerate {

inline size_t pulscan_truncate_fft_bins(size_t num_complex_bins) {
  return (num_complex_bins / 4096) * 4096;
}

inline size_t pulscan_select_samples_per_series(size_t available_samples,
                                                size_t planned_samples) {
  if(planned_samples == 0 || planned_samples > available_samples) {
    return available_samples;
  }
  return planned_samples;
}

inline bool pulscan_ensure_directory_exists(const std::string& path) {
  if(path.empty() || path == ".") {
    return true;
  }

  std::string current;
  size_t position = 0;
  if(path[0] == '/') {
    current = "/";
    position = 1;
  }

  while(position <= path.size()) {
    const size_t next_separator = path.find('/', position);
    const std::string part = path.substr(position, next_separator - position);
    if(!part.empty() && part != ".") {
      if(!current.empty() && current.back() != '/') {
        current += "/";
      }
      current += part;

      struct stat path_stat;
      if(stat(current.c_str(), &path_stat) != 0) {
        if(mkdir(current.c_str(), 0755) != 0 && errno != EEXIST) {
          return false;
        }
      }
      else if(!S_ISDIR(path_stat.st_mode)) {
        return false;
      }
    }

    if(next_separator == std::string::npos) {
      break;
    }
    position = next_separator + 1;
  }

  return true;
}

inline std::string pulscan_join_path(const std::string& directory,
                                     const std::string& filename) {
  if(directory.empty() || directory == ".") {
    return filename;
  }

  if(directory.back() == '/') {
    return directory + filename;
  }

  return directory + "/" + filename;
}

inline std::string pulscan_fft_dump_directory(const char* dump_setting) {
  if(!dump_setting || dump_setting[0] == '\0' ||
     std::strcmp(dump_setting, "1") == 0) {
    return "pulscan_fft_dumps";
  }

  return std::string(dump_setting);
}

struct pulscan_workspace {
  float*             d_time_series;
  size_t             d_time_series_capacity;
  float2*            d_fft_output;
  size_t             d_fft_capacity;
  float*             d_real;
  size_t             d_real_capacity;
  float*             d_imag;
  size_t             d_imag_capacity;
  float*             d_power;
  size_t             d_power_capacity;
  float*             d_decimated2;
  size_t             d_decimated2_capacity;
  float*             d_decimated3;
  size_t             d_decimated3_capacity;
  float*             d_decimated4;
  size_t             d_decimated4_capacity;
  pulscan_candidate* d_candidates_h1;
  size_t             d_candidates_h1_capacity;
  pulscan_candidate* d_candidates_h2;
  size_t             d_candidates_h2_capacity;
  pulscan_candidate* d_candidates_h3;
  size_t             d_candidates_h3_capacity;
  pulscan_candidate* d_candidates_h4;
  size_t             d_candidates_h4_capacity;
  pulscan_candidate* h_candidate_buffer;
  size_t             h_candidate_capacity;
  float*             h_fft_real;
  float*             h_fft_imag;
  size_t             h_fft_capacity;
  cufftHandle        fft_plan;
  size_t             fft_length;
  int                fft_batch;
  bool               fft_plan_valid;
  float*             h_pinned_batch;
  size_t             h_pinned_capacity;

  pulscan_workspace() :
      d_time_series(nullptr),
      d_time_series_capacity(0),
      d_fft_output(nullptr),
      d_fft_capacity(0),
      d_real(nullptr),
      d_real_capacity(0),
      d_imag(nullptr),
      d_imag_capacity(0),
      d_power(nullptr),
      d_power_capacity(0),
      d_decimated2(nullptr),
      d_decimated2_capacity(0),
      d_decimated3(nullptr),
      d_decimated3_capacity(0),
      d_decimated4(nullptr),
      d_decimated4_capacity(0),
      d_candidates_h1(nullptr),
      d_candidates_h1_capacity(0),
      d_candidates_h2(nullptr),
      d_candidates_h2_capacity(0),
      d_candidates_h3(nullptr),
      d_candidates_h3_capacity(0),
      d_candidates_h4(nullptr),
      d_candidates_h4_capacity(0),
      h_candidate_buffer(nullptr),
      h_candidate_capacity(0),
      h_fft_real(nullptr),
      h_fft_imag(nullptr),
      h_fft_capacity(0),
      fft_plan(0),
      fft_length(0),
      fft_batch(0),
      fft_plan_valid(false),
      h_pinned_batch(nullptr),
      h_pinned_capacity(0) {}

  ~pulscan_workspace() { release(); }

  bool ensure_device(size_t elements) {
    if(elements <= d_time_series_capacity) {
      return true;
    }
    if(d_time_series) {
      cudaFree(d_time_series);
      d_time_series = nullptr;
      d_time_series_capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&d_time_series),
                                 elements * sizeof(float));
    if(err != cudaSuccess) {
      return false;
    }
    d_time_series_capacity = elements;
    return true;
  }

  bool ensure_fft_output(size_t elements) {
    if(elements <= d_fft_capacity) {
      return true;
    }
    if(d_fft_output) {
      cudaFree(d_fft_output);
      d_fft_output = nullptr;
      d_fft_capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&d_fft_output),
                                 elements * sizeof(float2));
    if(err != cudaSuccess) {
      return false;
    }
    d_fft_capacity = elements;
    return true;
  }

  cufftResult ensure_fft_plan(int length, int batch) {
    if(fft_plan_valid && static_cast<size_t>(length) == fft_length &&
       batch == fft_batch) {
      return CUFFT_SUCCESS;
    }
    if(fft_plan_valid) {
      cufftDestroy(fft_plan);
      fft_plan_valid = false;
    }
    cufftResult plan_status = cufftPlan1d(&fft_plan, length, CUFFT_R2C, batch);
    if(plan_status != CUFFT_SUCCESS) {
      fft_plan = 0;
      fft_length = 0;
      fft_batch = 0;
      return plan_status;
    }
    fft_length = static_cast<size_t>(length);
    fft_batch = batch;
    fft_plan_valid = true;
    return CUFFT_SUCCESS;
  }

  bool ensure_real(size_t elements) {
    if(elements <= d_real_capacity) {
      return true;
    }
    if(d_real) {
      cudaFree(d_real);
      d_real = nullptr;
      d_real_capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&d_real),
                                 elements * sizeof(float));
    if(err != cudaSuccess) {
      return false;
    }
    d_real_capacity = elements;
    return true;
  }

  bool ensure_imag(size_t elements) {
    if(elements <= d_imag_capacity) {
      return true;
    }
    if(d_imag) {
      cudaFree(d_imag);
      d_imag = nullptr;
      d_imag_capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&d_imag),
                                 elements * sizeof(float));
    if(err != cudaSuccess) {
      return false;
    }
    d_imag_capacity = elements;
    return true;
  }

  bool ensure_power(size_t elements) {
    if(elements <= d_power_capacity) {
      return true;
    }
    if(d_power) {
      cudaFree(d_power);
      d_power = nullptr;
      d_power_capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&d_power),
                                 elements * sizeof(float));
    if(err != cudaSuccess) {
      return false;
    }
    d_power_capacity = elements;
    return true;
  }

  bool ensure_decimated_buffer(float*& buffer,
                               size_t& capacity,
                               size_t elements) {
    if(elements <= capacity) {
      return true;
    }
    if(buffer) {
      cudaFree(buffer);
      buffer = nullptr;
      capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&buffer),
                                 elements * sizeof(float));
    if(err != cudaSuccess) {
      return false;
    }
    capacity = elements;
    return true;
  }

  bool ensure_candidates(size_t elements) {
    if(elements <= d_candidates_h1_capacity) {
      return true;
    }
    if(d_candidates_h1) {
      cudaFree(d_candidates_h1);
      d_candidates_h1 = nullptr;
      d_candidates_h1_capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&d_candidates_h1),
                                 elements * sizeof(pulscan_candidate));
    if(err != cudaSuccess) {
      return false;
    }
    d_candidates_h1_capacity = elements;
    return true;
  }

  bool ensure_candidates_harmonic(pulscan_candidate*& buffer,
                                  size_t& capacity,
                                  size_t elements) {
    if(elements <= capacity) {
      return true;
    }
    if(buffer) {
      cudaFree(buffer);
      buffer = nullptr;
      capacity = 0;
    }
    cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&buffer),
                                 elements * sizeof(pulscan_candidate));
    if(err != cudaSuccess) {
      return false;
    }
    capacity = elements;
    return true;
  }

  bool ensure_host(size_t elements) {
    if(elements <= h_pinned_capacity) {
      return true;
    }
    if(h_pinned_batch) {
      cudaFreeHost(h_pinned_batch);
      h_pinned_batch = nullptr;
      h_pinned_capacity = 0;
    }
    cudaError_t err = cudaMallocHost(reinterpret_cast<void**>(&h_pinned_batch),
                                     elements * sizeof(float));
    if(err != cudaSuccess) {
      return false;
    }
    h_pinned_capacity = elements;
    return true;
  }

  bool ensure_host_candidates(size_t elements) {
    if(elements <= h_candidate_capacity) {
      return true;
    }
    if(h_candidate_buffer) {
      free(h_candidate_buffer);
      h_candidate_buffer = nullptr;
      h_candidate_capacity = 0;
    }
    h_candidate_buffer = static_cast<pulscan_candidate*>(
        malloc(elements * sizeof(pulscan_candidate)));
    if(!h_candidate_buffer) {
      return false;
    }
    h_candidate_capacity = elements;
    return true;
  }

  bool ensure_host_fft(size_t bytes) {
    size_t elements = bytes / sizeof(float);
    if(elements <= h_fft_capacity) {
      return true;
    }
    free(h_fft_real);
    free(h_fft_imag);
    h_fft_real = static_cast<float*>(malloc(elements * sizeof(float)));
    h_fft_imag = static_cast<float*>(malloc(elements * sizeof(float)));
    if(!h_fft_real || !h_fft_imag) {
      free(h_fft_real);
      free(h_fft_imag);
      h_fft_real = nullptr;
      h_fft_imag = nullptr;
      h_fft_capacity = 0;
      return false;
    }
    h_fft_capacity = elements;
    return true;
  }

  void release() {
    if(d_time_series) {
      cudaFree(d_time_series);
      d_time_series = nullptr;
    }
    d_time_series_capacity = 0;
    if(d_fft_output) {
      cudaFree(d_fft_output);
      d_fft_output = nullptr;
    }
    d_fft_capacity = 0;
    if(d_real) {
      cudaFree(d_real);
      d_real = nullptr;
    }
    d_real_capacity = 0;
    if(d_imag) {
      cudaFree(d_imag);
      d_imag = nullptr;
    }
    d_imag_capacity = 0;
    if(d_power) {
      cudaFree(d_power);
      d_power = nullptr;
    }
    d_power_capacity = 0;
    if(d_decimated2) {
      cudaFree(d_decimated2);
      d_decimated2 = nullptr;
    }
    d_decimated2_capacity = 0;
    if(d_decimated3) {
      cudaFree(d_decimated3);
      d_decimated3 = nullptr;
    }
    d_decimated3_capacity = 0;
    if(d_decimated4) {
      cudaFree(d_decimated4);
      d_decimated4 = nullptr;
    }
    d_decimated4_capacity = 0;
    if(d_candidates_h1) {
      cudaFree(d_candidates_h1);
      d_candidates_h1 = nullptr;
    }
    d_candidates_h1_capacity = 0;
    if(d_candidates_h2) {
      cudaFree(d_candidates_h2);
      d_candidates_h2 = nullptr;
    }
    d_candidates_h2_capacity = 0;
    if(d_candidates_h3) {
      cudaFree(d_candidates_h3);
      d_candidates_h3 = nullptr;
    }
    d_candidates_h3_capacity = 0;
    if(d_candidates_h4) {
      cudaFree(d_candidates_h4);
      d_candidates_h4 = nullptr;
    }
    d_candidates_h4_capacity = 0;
    if(fft_plan_valid) {
      cufftDestroy(fft_plan);
      fft_plan_valid = false;
    }
    fft_length = 0;
    fft_batch = 0;
    if(h_pinned_batch) {
      cudaFreeHost(h_pinned_batch);
      h_pinned_batch = nullptr;
    }
    h_pinned_capacity = 0;
    if(h_candidate_buffer) {
      free(h_candidate_buffer);
      h_candidate_buffer = nullptr;
    }
    h_candidate_capacity = 0;
    if(h_fft_real) {
      free(h_fft_real);
      h_fft_real = nullptr;
    }
    if(h_fft_imag) {
      free(h_fft_imag);
      h_fft_imag = nullptr;
    }
    h_fft_capacity = 0;
  }
};

class aa_pulscan_runner {
public:
  const std::vector<pulscan_host_candidate>& candidates() const {
    return m_candidates;
  }

  void clear() { m_candidates.clear(); }

  void release() {
    m_candidates.clear();
    m_workspace.release();
  }

  bool run(aa_periodicity_strategy& periodicity_strategy,
           float***                 output_buffer,
           const std::vector<int>&  in_bin,
           long                     inc,
           int&                     pipeline_error,
           TimeLog&                 time_log) {
    if(pipeline_error != PIPELINE_ERROR_NO_ERROR) {
      return false;
    }

    aa_gpu_timer timer;
    timer.Start();
    m_candidates.clear();

    struct pulscan_batch_plan {
      int    range_index;
      int    batch_index;
      int    dm_offset;
      int    dm_count;
      size_t samples_per_series;
      double dm_low;
      double dm_step;
      double sampling_time;
      int    in_bin;
    };

    std::vector<pulscan_batch_plan> batch_plan;
    batch_plan.reserve(
        static_cast<size_t>(std::max(0, periodicity_strategy.nRanges())));

    size_t total_dm_series = 0;
    size_t min_samples_per_series = std::numeric_limits<size_t>::max();
    size_t max_samples_per_series = 0;

    const int n_periodicity_ranges = periodicity_strategy.nRanges();
    for(int r = 0; r < n_periodicity_ranges; ++r) {
      aa_periodicity_range current_range =
          periodicity_strategy.get_periodicity_range(static_cast<size_t>(r));
      const int in_bin_factor = std::max(1, in_bin[static_cast<size_t>(r)]);
      const size_t samples_per_dm = static_cast<size_t>(inc / in_bin_factor);
      if(samples_per_dm == 0) {
        continue;
      }

      const size_t batch_count = current_range.batches.size();
      for(size_t b = 0; b < batch_count; ++b) {
        const aa_periodicity_batch& batch = current_range.batches[b];
        const int dm_count = static_cast<int>(batch.nDMs_per_batch);
        if(dm_count <= 0) {
          continue;
        }

        // The periodicity batch size is budgeted against this corrected
        // length, not the full dedispersed series length.
        const size_t samples_for_series =
            pulscan_select_samples_per_series(samples_per_dm,
                                              batch.nTimesamples_to_copy);

        pulscan_batch_plan plan_entry{};
        plan_entry.range_index = r;
        plan_entry.batch_index = static_cast<int>(b);
        plan_entry.dm_offset = static_cast<int>(batch.DM_shift);
        plan_entry.dm_count = dm_count;
        plan_entry.samples_per_series = samples_for_series;
        plan_entry.dm_low = current_range.range.dm_low();
        plan_entry.dm_step = current_range.range.dm_step();
        plan_entry.sampling_time = current_range.range.sampling_time();
        plan_entry.in_bin = in_bin_factor;
        batch_plan.push_back(plan_entry);

        total_dm_series += static_cast<size_t>(dm_count);
        min_samples_per_series =
            std::min(min_samples_per_series, samples_for_series);
        max_samples_per_series =
            std::max(max_samples_per_series, samples_for_series);
      }
    }

    std::printf("Pulscan staging prepared %zu batches and %zu DM series.\n",
                batch_plan.size(),
                total_dm_series);
    if(!batch_plan.empty()) {
      std::printf(
          "Pulscan series length range: [%zu, %zu] samples (bin-corrected).\n",
          min_samples_per_series,
          max_samples_per_series);
    }
    else {
      std::printf("Pulscan staging skipped: no qualified DM series found.\n");
    }

    cudaError_t cuda_status = cudaSuccess;
    for(const auto& plan_entry : batch_plan) {
      if(pipeline_error != PIPELINE_ERROR_NO_ERROR) {
        break;
      }

      const size_t dm_series =
          static_cast<size_t>(std::max(0, plan_entry.dm_count));
      const size_t samples_per_series = plan_entry.samples_per_series;
      const size_t batch_elements = dm_series * samples_per_series;
      if(batch_elements == 0) {
        continue;
      }

      if(!m_workspace.ensure_host(batch_elements) ||
         !m_workspace.ensure_device(batch_elements)) {
        LOG(log_level::error, "Pulscan workspace allocation failed.");
        pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
        break;
      }

      float* host_buffer = m_workspace.h_pinned_batch;
      const size_t host_stride = samples_per_series;
      for(int dm_idx = 0; dm_idx < plan_entry.dm_count; ++dm_idx) {
        float* series_ptr =
            output_buffer[plan_entry.range_index][plan_entry.dm_offset + dm_idx];
        if(series_ptr == nullptr) {
          continue;
        }
        std::copy_n(series_ptr,
                    host_stride,
                    host_buffer + static_cast<size_t>(dm_idx) * host_stride);
      }

      cuda_status = cudaMemcpy(m_workspace.d_time_series,
                               host_buffer,
                               batch_elements * sizeof(float),
                               cudaMemcpyHostToDevice);
      if(cuda_status != cudaSuccess) {
        LOG(log_level::error,
            "Pulscan: failed to copy time series batch to device (" +
                std::string(cudaGetErrorString(cuda_status)) + ")");
        pipeline_error = PIPELINE_ERROR_COPY_TO_DEVICE;
        break;
      }

      std::printf(
          "Pulscan staged range %d batch %d: %d DM trials, %zu samples each.\n",
          plan_entry.range_index,
          plan_entry.batch_index,
          plan_entry.dm_count,
          host_stride);

      const int fft_length = static_cast<int>(samples_per_series);
      const int fft_batch = static_cast<int>(dm_series);
      if(fft_length <= 1 || fft_batch <= 0) {
        continue;
      }

      const cufftResult plan_status =
          m_workspace.ensure_fft_plan(fft_length, fft_batch);
      if(plan_status != CUFFT_SUCCESS) {
        LOG(log_level::error,
            "Pulscan: failed to create cuFFT plan (length=" +
                std::to_string(fft_length) + ", batch=" +
                std::to_string(fft_batch) + ", status=" +
                std::to_string(static_cast<int>(plan_status)) + ")");
        pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
        break;
      }

      const size_t freq_bins = static_cast<size_t>(fft_length / 2 + 1);
      const size_t fft_elements = dm_series * freq_bins;
      if(!m_workspace.ensure_fft_output(fft_elements)) {
        LOG(log_level::error,
            "Pulscan: failed to allocate FFT output buffer for " +
                std::to_string(fft_elements) + " complex samples.");
        pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
        break;
      }

      cufftResult fft_status = cufftExecR2C(
          m_workspace.fft_plan,
          reinterpret_cast<cufftReal*>(m_workspace.d_time_series),
          reinterpret_cast<cufftComplex*>(m_workspace.d_fft_output));
      if(fft_status != CUFFT_SUCCESS) {
        LOG(log_level::error,
            "Pulscan: cuFFT execution failed (status=" +
                std::to_string(fft_status) + ")");
        pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
        break;
      }

      std::printf(
          "Pulscan FFT completed for range %d batch %d (%d series, %zu frequency bins).\n",
          plan_entry.range_index,
          plan_entry.batch_index,
          plan_entry.dm_count,
          freq_bins);

      const size_t usable_freq_bins = pulscan_truncate_fft_bins(freq_bins);
      if(usable_freq_bins == 0) {
        std::printf("Pulscan skipped range %d batch %d: %zu FFT bins truncate to zero.\n",
                    plan_entry.range_index,
                    plan_entry.batch_index,
                    freq_bins);
        continue;
      }

      if(!m_workspace.ensure_real(usable_freq_bins) ||
         !m_workspace.ensure_imag(usable_freq_bins) ||
         !m_workspace.ensure_power(usable_freq_bins)) {
        LOG(log_level::error,
            "Pulscan: failed to allocate magnitude workspace of size " +
                std::to_string(usable_freq_bins));
        pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
        break;
      }

      const float logp_threshold = -10.0f;
      const int boxcar_thread_count = 256;
      const int num_candidates_per_block = []() {
        int zmax = 256;
        int per_block = 1;
        int value = 1;
        while(value <= zmax) {
          value *= 2;
          per_block += 1;
        }
        return per_block;
      }();

      for(int dm_idx = 0; dm_idx < plan_entry.dm_count; ++dm_idx) {
        cudaStream_t stream = nullptr;
        float2* dm_fft =
            m_workspace.d_fft_output + static_cast<size_t>(dm_idx) * freq_bins;
        const dim3 fft_block(boxcar_thread_count);
        const dim3 fft_grid(
            (static_cast<int>(usable_freq_bins) + boxcar_thread_count - 1) /
            boxcar_thread_count);

        call_kernel_pulscan_separate_components(fft_grid,
                                                fft_block,
                                                stream,
                                                dm_fft,
                                                m_workspace.d_real,
                                                m_workspace.d_imag,
                                                static_cast<long>(usable_freq_bins));

        if(const char* dump_setting = std::getenv("AA_PULSCAN_DUMP_FFT")) {
          const size_t real_bytes = usable_freq_bins * sizeof(float);
          if(m_workspace.ensure_host_fft(real_bytes)) {
            cudaMemcpy(m_workspace.h_fft_real,
                       m_workspace.d_real,
                       real_bytes,
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(m_workspace.h_fft_imag,
                       m_workspace.d_imag,
                       real_bytes,
                       cudaMemcpyDeviceToHost);
            const std::string dump_dir =
                pulscan_fft_dump_directory(dump_setting);
            if(!pulscan_ensure_directory_exists(dump_dir)) {
              LOG(log_level::error,
                  "Pulscan: failed to create FFT dump directory " + dump_dir);
            }
            else {
              char filename[256];
              std::snprintf(filename,
                            sizeof(filename),
                            "pulscan_fft_dm%03d_part%d.bin",
                            plan_entry.range_index,
                            dm_idx);
              const std::string filepath =
                  pulscan_join_path(dump_dir, filename);
              FILE* fp_fft = std::fopen(filepath.c_str(), "wb");
              if(fp_fft) {
                for(size_t i = 0; i < usable_freq_bins; ++i) {
                  std::fwrite(&m_workspace.h_fft_real[i], sizeof(float), 1, fp_fft);
                  std::fwrite(&m_workspace.h_fft_imag[i], sizeof(float), 1, fp_fft);
                }
                std::fclose(fp_fft);
              }
            }
          }
        }

        const dim3 normalise_grid(
            static_cast<unsigned int>(usable_freq_bins / 4096));
        const dim3 normalise_block(1024);
        call_kernel_pulscan_median_normalisation(normalise_grid,
                                                 normalise_block,
                                                 stream,
                                                 m_workspace.d_real);
        call_kernel_pulscan_median_normalisation(normalise_grid,
                                                 normalise_block,
                                                 stream,
                                                 m_workspace.d_imag);

        call_kernel_pulscan_magnitude_squared(fft_grid,
                                              fft_block,
                                              stream,
                                              m_workspace.d_real,
                                              m_workspace.d_imag,
                                              m_workspace.d_power,
                                              static_cast<long>(usable_freq_bins));

        cuda_status = cudaGetLastError();
        if(cuda_status != cudaSuccess) {
          LOG(log_level::error,
              "Pulscan: magnitude kernels failed (" +
                  std::string(cudaGetErrorString(cuda_status)) + ")");
          pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
          break;
        }

        const size_t dec2_bins = usable_freq_bins / 2;
        const size_t dec3_bins = usable_freq_bins / 3;
        const size_t dec4_bins = usable_freq_bins / 4;
        if((dec2_bins > 0 &&
            !m_workspace.ensure_decimated_buffer(m_workspace.d_decimated2,
                                                 m_workspace.d_decimated2_capacity,
                                                 dec2_bins)) ||
           (dec3_bins > 0 &&
            !m_workspace.ensure_decimated_buffer(m_workspace.d_decimated3,
                                                 m_workspace.d_decimated3_capacity,
                                                 dec3_bins)) ||
           (dec4_bins > 0 &&
            !m_workspace.ensure_decimated_buffer(m_workspace.d_decimated4,
                                                 m_workspace.d_decimated4_capacity,
                                                 dec4_bins))) {
          LOG(log_level::error,
              "Pulscan: failed to allocate decimated power buffers.");
          pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
          break;
        }

        if(dec2_bins > 0 || dec3_bins > 0 || dec4_bins > 0) {
          const int decimation_grid = static_cast<int>(
              (static_cast<long>(dec2_bins) + boxcar_thread_count - 1) /
              boxcar_thread_count);
          if(decimation_grid > 0) {
            call_kernel_pulscan_decimate_harmonics(
                dim3(decimation_grid),
                dim3(boxcar_thread_count),
                stream,
                m_workspace.d_power,
                m_workspace.d_decimated2,
                m_workspace.d_decimated3,
                m_workspace.d_decimated4,
                static_cast<long>(usable_freq_bins));

            cuda_status = cudaGetLastError();
            if(cuda_status != cudaSuccess) {
              LOG(log_level::error,
                  "Pulscan: decimation kernel failed (" +
                      std::string(cudaGetErrorString(cuda_status)) + ")");
              pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
              break;
            }
          }
        }

        struct harmonic_run {
          int                num_harmonics;
          float*             device_input;
          size_t             input_length;
          pulscan_candidate* device_candidates;
          size_t             candidate_count;
          int                grid_size;
          int                num_sum;
          size_t             host_offset;
        };

        std::array<harmonic_run, 4> harmonic_runs{};
        size_t harmonic_run_count = 0;

        auto try_register_harmonic = [&](int                 num_harmonics,
                                         float*              input,
                                         size_t              length,
                                         pulscan_candidate*& buffer,
                                         size_t&             capacity) {
          if(pipeline_error != PIPELINE_ERROR_NO_ERROR || length == 0) {
            return;
          }
          const int grid_size = static_cast<int>(
              (static_cast<long>(length) + boxcar_thread_count - 1) /
              boxcar_thread_count);
          if(grid_size <= 0) {
            return;
          }
          const size_t candidate_count = static_cast<size_t>(grid_size) *
                                         static_cast<size_t>(num_candidates_per_block);
          if(candidate_count == 0) {
            return;
          }
          const bool ensured =
              (num_harmonics == 1)
                  ? m_workspace.ensure_candidates(candidate_count)
                  : m_workspace.ensure_candidates_harmonic(
                        buffer, capacity, candidate_count);
          if(!ensured) {
            LOG(log_level::error,
                "Pulscan: failed to allocate candidate buffer for harmonic " +
                    std::to_string(num_harmonics));
            pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
            return;
          }

          cuda_status =
              cudaMemset(buffer, 0, candidate_count * sizeof(pulscan_candidate));
          if(cuda_status != cudaSuccess) {
            LOG(log_level::error,
                "Pulscan: failed to zero candidate buffer (" +
                    std::string(cudaGetErrorString(cuda_status)) + ")");
            pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
            return;
          }

          harmonic_run& run = harmonic_runs[harmonic_run_count++];
          run.num_harmonics = num_harmonics;
          run.device_input = input;
          run.input_length = length;
          run.device_candidates = buffer;
          run.candidate_count = candidate_count;
          run.grid_size = grid_size;
          run.num_sum = (num_harmonics * (num_harmonics + 1)) / 2;
          run.host_offset = 0;
        };

        try_register_harmonic(1,
                              m_workspace.d_power,
                              usable_freq_bins,
                              m_workspace.d_candidates_h1,
                              m_workspace.d_candidates_h1_capacity);
        try_register_harmonic(2,
                              m_workspace.d_decimated2,
                              dec2_bins,
                              m_workspace.d_candidates_h2,
                              m_workspace.d_candidates_h2_capacity);
        try_register_harmonic(3,
                              m_workspace.d_decimated3,
                              dec3_bins,
                              m_workspace.d_candidates_h3,
                              m_workspace.d_candidates_h3_capacity);
        try_register_harmonic(4,
                              m_workspace.d_decimated4,
                              dec4_bins,
                              m_workspace.d_candidates_h4,
                              m_workspace.d_candidates_h4_capacity);

        if(pipeline_error != PIPELINE_ERROR_NO_ERROR) {
          break;
        }
        if(harmonic_run_count == 0) {
          continue;
        }

        for(size_t h = 0; h < harmonic_run_count; ++h) {
          harmonic_run& run = harmonic_runs[h];
          call_kernel_pulscan_boxcar_filter(dim3(run.grid_size),
                                            dim3(boxcar_thread_count),
                                            stream,
                                            run.device_input,
                                            run.device_candidates,
                                            run.num_harmonics,
                                            static_cast<long>(run.input_length),
                                            num_candidates_per_block);
        }

        for(size_t h = 0; h < harmonic_run_count; ++h) {
          harmonic_run& run = harmonic_runs[h];
          const int logp_grid = static_cast<int>(
              (static_cast<long>(run.candidate_count) + boxcar_thread_count - 1) /
              boxcar_thread_count);
          call_kernel_pulscan_calculate_logp(dim3(logp_grid),
                                             dim3(boxcar_thread_count),
                                             stream,
                                             run.device_candidates,
                                             static_cast<long>(run.candidate_count),
                                             run.num_sum);
        }

        cuda_status = cudaDeviceSynchronize();
        if(cuda_status != cudaSuccess) {
          LOG(log_level::error,
              "Pulscan: candidate kernels failed (" +
                  std::string(cudaGetErrorString(cuda_status)) + ")");
          pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
          break;
        }

        size_t total_candidates = 0;
        for(size_t h = 0; h < harmonic_run_count; ++h) {
          total_candidates += harmonic_runs[h].candidate_count;
        }
        if(total_candidates == 0) {
          continue;
        }

        if(!m_workspace.ensure_host_candidates(total_candidates)) {
          LOG(log_level::error,
              "Pulscan: failed to allocate host candidate buffer for " +
                  std::to_string(total_candidates) + " entries.");
          pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
          break;
        }

        size_t host_offset = 0;
        for(size_t h = 0; h < harmonic_run_count; ++h) {
          harmonic_run& run = harmonic_runs[h];
          run.host_offset = host_offset;
          if(run.candidate_count > 0) {
            cuda_status = cudaMemcpy(m_workspace.h_candidate_buffer + host_offset,
                                     run.device_candidates,
                                     run.candidate_count * sizeof(pulscan_candidate),
                                     cudaMemcpyDeviceToHost);
            if(cuda_status != cudaSuccess) {
              LOG(log_level::error,
                  "Pulscan: failed to copy candidates to host (" +
                      std::string(cudaGetErrorString(cuda_status)) + ")");
              pipeline_error = PIPELINE_ERROR_COPY_TO_HOST;
              break;
            }
          }
          host_offset += run.candidate_count;
        }
        if(pipeline_error != PIPELINE_ERROR_NO_ERROR) {
          break;
        }

        const double dm_value =
            plan_entry.dm_low +
            static_cast<double>(plan_entry.dm_offset + dm_idx) * plan_entry.dm_step;
        size_t candidates_added = 0;
        const double dt =
            plan_entry.sampling_time * static_cast<double>(plan_entry.in_bin);
        const double denom = static_cast<double>(fft_length) * dt;

        for(size_t h = 0; h < harmonic_run_count; ++h) {
          harmonic_run& run = harmonic_runs[h];
          pulscan_candidate* host_slice =
              m_workspace.h_candidate_buffer + run.host_offset;
          for(size_t c = 0; c < run.candidate_count; ++c) {
            const pulscan_candidate candidate = host_slice[c];
            if(candidate.logp >= logp_threshold || candidate.r == 0 ||
               candidate.z == 0) {
              continue;
            }
            const double frequency_hz =
                (denom > 0.0) ? static_cast<double>(candidate.r) / denom : 0.0;
            const double period_s = (frequency_hz > 0.0) ? 1.0 / frequency_hz : 0.0;

            pulscan_host_candidate host_candidate{};
            host_candidate.dm = dm_value;
            host_candidate.frequency_hz = frequency_hz;
            host_candidate.period_s = period_s;
            host_candidate.power = candidate.power;
            host_candidate.logp = candidate.logp;
            host_candidate.sigma = pulscan_logp_to_sigma(candidate.logp);
            host_candidate.r = candidate.r;
            host_candidate.z = candidate.z;
            host_candidate.numharm = candidate.numharm;
            host_candidate.range_index = plan_entry.range_index;
            host_candidate.dm_index = plan_entry.dm_offset + dm_idx;
            host_candidate.bin_index = candidate.r;
            if(!std::isfinite(host_candidate.sigma)) {
              host_candidate.sigma = 0.0f;
            }
            m_candidates.push_back(host_candidate);
            ++candidates_added;
          }
        }

        std::printf("Pulscan extracted %zu candidates for DM %.3f (dt=%f s).\n",
                    candidates_added,
                    dm_value,
                    dt);
      }

      if(pipeline_error != PIPELINE_ERROR_NO_ERROR) {
        break;
      }
    }

    std::sort(m_candidates.begin(),
              m_candidates.end(),
              [](const pulscan_host_candidate& a,
                 const pulscan_host_candidate& b) { return a.logp < b.logp; });

    std::printf("Pulscan accumulated %zu candidates so far.\n",
                m_candidates.size());

    timer.Stop();
    time_log.adding("Pulscan", "total", timer.Elapsed());
    time_log.adding("Total", "total", timer.Elapsed());
    return (pipeline_error == PIPELINE_ERROR_NO_ERROR);
  }

private:
  pulscan_workspace                   m_workspace;
  std::vector<pulscan_host_candidate> m_candidates;
};

} // namespace astroaccelerate

#endif // AA_ENABLE_PULSCAN

#endif // ASTRO_ACCELERATE_AA_PULSCAN_RUNNER_HPP
