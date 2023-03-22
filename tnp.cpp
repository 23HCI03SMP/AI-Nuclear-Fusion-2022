#include "include/traj.h"
void tnp(float *Ea, float *Ba,
         float *pos0x, float *pos0y, float *pos0z,
         float *pos1x, float *pos1y, float *pos1z,
         float cf[2],
         unsigned int ci[2])
{
   unsigned int n = n_partd;//ci[0];
   unsigned int nc = n_cells * ncoeff * 3; // trilin constatnts have 8 coefficients 3 components

   static bool fastIO;
   static bool first = true;
   if(first){ // get whether or not we are on an iGPU/similar, and can use certain memmory optimizations
      bool temp;
      default_device_g.getInfo(CL_DEVICE_HOST_UNIFIED_MEMORY, &temp);
      cout << "Using unified memory: " << temp << " ";
      fastIO = temp;
      first = false;
   }

   // create buffers on the device
   /** IMPORTANT: do not use CL_MEM_USE_HOST_PTR if on dGPU **/
   /** HOST_PTR is only used so that memory is not copied, but instead shared between CPU and iGPU in RAM**/
   // Note that special alignment has been given to Ea, Ba, y0, z0, x0, x1, y1 in order to actually do this properly

   // Assume buffers A, B, I, J (Ea, Ba, ci, cf) will always be the same. Then we save a bit of time.
   static cl::Buffer buffer_A(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_ONLY, sizeof(float) * nc, fastIO ? Ea : NULL);
   static cl::Buffer buffer_B(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_ONLY, sizeof(float) * nc, fastIO ? Ba : NULL);
   cl::Buffer buffer_C(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_WRITE, sizeof(float) * n, fastIO ? pos0x : NULL);  // x0
   cl::Buffer buffer_D(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_WRITE, sizeof(float) * n, fastIO ? pos0y : NULL);  // y0
   cl::Buffer buffer_E(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_WRITE, sizeof(float) * n, fastIO ? pos0z : NULL);  // z0
   cl::Buffer buffer_F(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_WRITE, sizeof(float) * n, fastIO ? pos1x : NULL);  // x1
   cl::Buffer buffer_G(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_WRITE, sizeof(float) * n, fastIO ? pos1y : NULL);  // y1
   cl::Buffer buffer_H(context_g, (fastIO ? CL_MEM_USE_HOST_PTR : 0) | CL_MEM_READ_WRITE, sizeof(float) * n, fastIO ? pos1z : NULL);  // z1

   // create queue to which we will push commands for the device.
   static cl::CommandQueue queue(context_g, default_device_g);
   cl::Kernel kernel_add = cl::Kernel(program_g, "tnp_k_implicit"); // select the kernel program to run
   // write input arrays to the device
   if(fastIO){//is mapping required? // Yes we might need to map because OpenCL does not guarantee that the data will be shared, alternatively use SVM
   //auto * mapped_buffer_C = (float *)queue.enqueueMapBuffer(buffer_C, CL_TRUE, CL_MAP_WRITE, 0, sizeof(float) * n); queue.enqueueUnmapMemObject(buffer_C, mapped_buffer_C);
   } else {
      queue.enqueueWriteBuffer(buffer_A, CL_TRUE, 0, sizeof(float) * nc, Ea);
      queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(float) * nc, Ba);
      queue.enqueueWriteBuffer(buffer_C, CL_TRUE, 0, sizeof(float) * n, pos0x);
      queue.enqueueWriteBuffer(buffer_D, CL_TRUE, 0, sizeof(float) * n, pos0y);
      queue.enqueueWriteBuffer(buffer_E, CL_TRUE, 0, sizeof(float) * n, pos0z);
      queue.enqueueWriteBuffer(buffer_F, CL_TRUE, 0, sizeof(float) * n, pos1x);
      queue.enqueueWriteBuffer(buffer_G, CL_TRUE, 0, sizeof(float) * n, pos1y);
      queue.enqueueWriteBuffer(buffer_H, CL_TRUE, 0, sizeof(float) * n, pos1z);   
   }
   //return;
   // set arguments to be fed into the kernel program
   kernel_add.setArg(0, buffer_A);                       // the 1st argument to the kernel program Ea
   kernel_add.setArg(1, buffer_B);                       // Ba
   kernel_add.setArg(2, buffer_C);                       // x0
   kernel_add.setArg(3, buffer_D);                       // y0
   kernel_add.setArg(4, buffer_E);                       // z0
   kernel_add.setArg(5, buffer_F);                       // x1
   kernel_add.setArg(6, buffer_G);                       // y1
   kernel_add.setArg(7, buffer_H);                       // z1
   kernel_add.setArg(8, sizeof(float), &cf[0]); // Bconst
   kernel_add.setArg(9, sizeof(float), &cf[1]); // Econst
   kernel_add.setArg(10, sizeof(int), &ci[0]); // npart
   kernel_add.setArg(11, sizeof(int), &ci[1]); // ncalc

   // run the kernel
   queue.enqueueNDRangeKernel(kernel_add, cl::NullRange, cl::NDRange(n), cl::NullRange);

   // read result arrays from the device to main memory
   if (fastIO){ //is mapping required?
   //mapped_buffer_C = (float *)queue.enqueueMapBuffer(buffer_C, CL_TRUE, CL_MAP_READ, 0, sizeof(float) * n); queue.enqueueUnmapMemObject(buffer_C, mapped_buffer_C);
   } else {
      queue.enqueueReadBuffer(buffer_C, CL_TRUE, 0, sizeof(float) * n, pos0x);
      queue.enqueueReadBuffer(buffer_D, CL_TRUE, 0, sizeof(float) * n, pos0y);
      queue.enqueueReadBuffer(buffer_E, CL_TRUE, 0, sizeof(float) * n, pos0z);
      queue.enqueueReadBuffer(buffer_F, CL_TRUE, 0, sizeof(float) * n, pos1x);
      queue.enqueueReadBuffer(buffer_G, CL_TRUE, 0, sizeof(float) * n, pos1y);
      queue.enqueueReadBuffer(buffer_H, CL_TRUE, 0, sizeof(float) * n, pos1z);
   }
   queue.finish(); // wait for the end of the kernel program
}