
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
//#include <helper_cuda.h>

#include <stdio.h>
#include <assert.h>

// The launch configurator returned block size
int blockSize;
// The actual grid size needed, based on input size
int gridSize;

cudaError_t multiplyWithCuda(const double* A, const double* B, double* C, int numElements);
cudaError_t divideWithCuda(const double* A, const double* B, double* C, int numElements);
cudaError_t addWithCuda(const double* A, const double* B, double* C, int numElements);
cudaError_t subtractWithCuda(const double* A, const double* B, double* C, int numElements);
cudaError_t MATmultiplyWithCuda(const double* A, const double* B, double* C, int numRows, int numCols);

__global__ void vectorMultiply(const double* A, const double* B, double* C, int numElements) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements) {
        C[i] = A[i] * B[i];
    }
}

__global__ void vectorDivide(const double* A, const double* B, double* C, int numElements) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;


    if (i < numElements) {
        C[i] = A[i] / B[i];
    }
}

__global__ void vectorAdd(const double* A, const double* B, double* C, int numElements) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements) {
        C[i] = A[i] + B[i];
    }
}

__global__ void vectorSubtract(const double* A, const double* B, double* C, int numElements) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements) {
        C[i] = A[i] - B[i];
    }
}

__global__ void matrixMultiply(const double* A, const double* B, double* C, int N, int M, int N2, int M2) {
    
    int ROW = blockIdx.x * blockDim.x + threadIdx.x;
    int COL = blockIdx.y * blockDim.y + threadIdx.y;
    
    //printf("%i %i,%i %i,%i %i,%i %i\n", gridDim.x, gridDim.y, blockIdx.x, blockIdx.y, blockDim.x, blockDim.y, threadIdx.x, threadIdx.y);
    printf("%i %i\n", ROW, COL);
    double tmpSum = 0;

    if (ROW >= N || COL >= M2) {
        //printf("hmmm\n");
        return;
    }
    // each thread computes one element of the block sub-matrix
    for (int i = 0; i < N; i++) {
        tmpSum += A[ROW * M + i] * B[i * M2 + COL];
        //printf("%i %i\n", A[ROW * M + i], B[i * M + COL]);
    }
    C[ROW * M2 + COL] = tmpSum;


}

int notmain()
{
    const int arraySize = 5;
    const double a[arraySize] = { 1, 2, 3, 4, 5 };
    const double b[arraySize] = { 10, 20, 30, 40, 50 };
    double c[arraySize] = { 0 };

    // Add vectors in parallel.
    cudaError_t cudaStatus = multiplyWithCuda(a, b, c, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "multiplyWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} * {10,20,30,40,50} = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);
    
    cudaStatus = divideWithCuda(a, b, c, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "divideWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} / {10,20,30,40,50} = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);

    cudaStatus = addWithCuda(a, b, c, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);

    cudaStatus = subtractWithCuda(a, b, c, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "subtractWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} - {10,20,30,40,50} = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to multiply vectors in parallel.
cudaError_t multiplyWithCuda(const double* A, const double* B, double* C, int numElements)
{
    double*dev_a = 0;
    double*dev_b = 0;
    double*dev_c = 0;
    cudaError_t cudaStatus;
    

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, A, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, B, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with multiple grids/blocks.
    vectorMultiply<<<gridSize, blockSize >>>(dev_a, dev_b, dev_c, numElements);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "multiplyWithCuda launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(C, dev_c, numElements * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);
    
    return cudaStatus;
}


// Helper function for using CUDA to multiply vectors in parallel.
cudaError_t divideWithCuda(const double* A, const double* B, double* C, int numElements) {
    double* dev_a = 0;
    double* dev_b = 0;
    double* dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, A, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, B, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    vectorDivide<<<1, numElements >>> (dev_a, dev_b, dev_c, numElements);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "divideWithCuda launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(C, dev_c, numElements * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);

    return cudaStatus;
}


// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(const double* A, const double* B, double* C, int numElements) {
    double* dev_a = 0;
    double* dev_b = 0;
    double* dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, A, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, B, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    vectorAdd<<<1, numElements >>> (dev_a, dev_b, dev_c, numElements);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(C, dev_c, numElements * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);

    return cudaStatus;
}


// Helper function for using CUDA to subtract vectors in parallel.
cudaError_t subtractWithCuda(const double* A, const double* B, double* C, int numElements) {
    double* dev_a = 0;
    double* dev_b = 0;
    double* dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, numElements * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, A, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, B, numElements * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    vectorSubtract<<<1, numElements >>> (dev_a, dev_b, dev_c, numElements);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "subtractWithCuda launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(C, dev_c, numElements * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);

    return cudaStatus;
}

void cudaCheck() {
    struct cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, 0);
    fprintf(stdout, "using %i multiprocessors\n", properties.multiProcessorCount);
    fprintf(stdout, "max threads per processor: %i\n", properties.maxThreadsPerMultiProcessor);
    fprintf(stdout, "number of concurrent jobs %i\n", properties.multiProcessorCount * properties.maxThreadsPerMultiProcessor);

    



}


// Helper function for using CUDA to multiply vectors in parallel.
cudaError_t MATmultiplyWithCuda(const double* A, const double* B, double* C, int numRows_a, int numCols_a, int numRows_b, int numCols_b) {
    double* dev_a = 0;
    double* dev_b = 0;
    double* dev_c = 0;
    cudaError_t cudaStatus;
    int numElements_a = numRows_a * numCols_a;
    int numElements_b = numRows_b * numCols_b;
    int numElements_c = numRows_a * numCols_b;
    assert(numRows_b == numCols_a);

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, numElements_c * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, numElements_a * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, numElements_b * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, A, numElements_a * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, B, numElements_b * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with multiple grids/blocks.
    matrixMultiply << <1, dim3(numRows_a, numCols_b) >> > (dev_a, dev_b, dev_c, numRows_a, numCols_a, numRows_b, numCols_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "multiplyWithCuda launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(C, dev_c, numElements_c * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);

    return cudaStatus;
}





extern "C" {
    void wrapper(double*a, double* b, double* c, int numElements) {
        cudaError_t cudaStatus = multiplyWithCuda(a, b, c, numElements);
        assert(cudaStatus == cudaSuccess);
    }
    void cudaStart(int numElements) {
        int minGridSize;    // The minimum grid size needed to achieve the
                            // maximum occupancy for a full device
                            // launch
        
        cudaOccupancyMaxPotentialBlockSize(
            &minGridSize,
            &blockSize,
            (void*)vectorMultiply,
            0,
            numElements);

        // Round up according to array size
        gridSize = (numElements + blockSize - 1) / blockSize;
        fprintf(stdout, "mesh %i, grid %i, block %i\n", numElements, gridSize, blockSize);
        

        /*
        int d = 3, e = 2;
        int f = 2, g = 4;
        int tot = d * e;
        const double a[6] = { 1,1,2,2,3,3 };
        const double b[8] = { 1,1,1,1,2,2,2,2 };
        double c[12] = { 0 };

        cudaError_t error = MATmultiplyWithCuda(a, b, c, d, e, f, g);
        
        assert(error == cudaSuccess);

        fprintf(stdout, "\n\n%f %f %f\n%f %f %f\n%f %f %f\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8]);
        */
    }
    void cudaEnd() {
        cudaError_t cudaStatus = cudaDeviceReset();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceReset failed!");
        }
        fprintf(stdout, "Cuda Freed\n");
    }
}