#include "include/traj.h"
#include <string>
#include <fstream>
#include <streambuf>
 cl::Context context_g;
 cl::Device default_device_g;
 cl::Program program_g;

stringstream cl_build_options;
void add_build_option(string name, string param){
    cl_build_options << "-D" << name << "=" << param << " ";
}
void add_build_option(string name, int param) { add_build_option(name, to_string(param)); }
void add_build_option(string name, float param) {
        cl_build_options << "-D" << name << "=" << param << "f ";
}

void cl_set_build_options(float posL[3], float posH[3], float dd[3]){
    add_build_option("XLOW", posL[0]);
    add_build_option("YLOW", posL[1]);
    add_build_option("ZLOW", posL[2]);
    add_build_option("XHIGH", posH[0]);
    add_build_option("YHIGH", posH[1]);
    add_build_option("ZHIGH", posH[2]);
    add_build_option("DX", dd[0]);
    add_build_option("DY", dd[1]);
    add_build_option("DZ", dd[2]);
    add_build_option("NX", n_space_divx);
    add_build_option("NY", n_space_divy);
    add_build_option("NZ", n_space_divz);
    //add_build_option("NC", n_cells);
}
//void cl_start(cl::Context &context1, cl::Device &default_device1, cl::Program &program1)
void cl_start()
{
    int AA[1] = {-1};
#pragma omp target
    {
        AA[0] = omp_is_initial_device();
    }
    if (!AA[0])
    {
        cout << "Able to use GPU offloading with OMP!\n";
    }
    else
    {
        cout << "\nNo GPU on OMP\n";
    }
    // get all platforms (drivers)
    cl::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    cl::vector<cl::Device> devices;
    int platform_id = 0;
    int device_id = 0;

    std::cout << "Number of Platforms: " << platforms.size() << std::endl;

    for (cl::vector<cl::Platform>::iterator it = platforms.begin(); it != platforms.end(); ++it)
    {
        cl::Platform platform(*it);

        std::cout << "Platform ID: " << platform_id++ << std::endl;
        std::cout << "Platform Name: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
        std::cout << "Platform Vendor: " << platform.getInfo<CL_PLATFORM_VENDOR>() << std::endl;

        //     platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
        platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
        //       platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
        for (cl::vector<cl::Device>::iterator it2 = devices.begin(); it2 != devices.end(); ++it2)
        {
            cl::Device device(*it2);
            std::cout << "Number of Devices: " << devices.size() << std::endl;
            std::cout << "\tDevice " << device_id++ << ": " << std::endl;
            std::cout << "\t\tDevice Name: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
            std::cout << "\t\tDevice Type: " << device.getInfo<CL_DEVICE_TYPE>();
            std::cout << " (GPU: " << CL_DEVICE_TYPE_GPU << ", CPU: " << CL_DEVICE_TYPE_CPU << ")" << std::endl;
            std::cout << "\t\tDevice Vendor: " << device.getInfo<CL_DEVICE_VENDOR>() << std::endl;
            std::cout << "\t\tDevice Max Compute Units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
            std::cout << "\t\tDevice Global Memory: MB " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / 1024 / 1024 << std::endl;
            std::cout << "\t\tDevice Max Clock Frequency: MHz " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
            std::cout << "\t\tDevice Max Allocateable Memory MB: " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() / 1024 / 1024 << std::endl;
            std::cout << "\t\tDevice Local Memory: kB " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() / 1024 << std::endl;
            std::cout << "\t\tDevice Available: " << device.getInfo<CL_DEVICE_AVAILABLE>() << std::endl;
        }
        std::cout << std::endl;
    }

    cl::Platform::get(&platforms);
    cl::Platform default_platform = platforms[0];
    std::cout << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
    cl::Device default_device = devices[0];
    std::cout << "\t\tDevice Name: " << default_device.getInfo<CL_DEVICE_NAME>() << "\n";
    std::cout << "OpenCL Version: " << default_device.getInfo<CL_DEVICE_VERSION>() << std::endl;

    cl::Context context({default_device});

    cl::Program::Sources sources;

    // read in kernel code which calculates the next position
    std::ifstream t("cl_kernel_code.cl");
    std::string kernel_code((std::istreambuf_iterator<char>(t)),
                            std::istreambuf_iterator<char>());
    //    kernel_code << t.rdbuf();
    // std::string kernel_code ="";

    sources.push_back({kernel_code.c_str(), kernel_code.length()});

    cl::Program program(context, sources);
    //cout << cl_build_options.str() << endl;
    //exit(0);
    cl_int cl_err=program.build({default_device}, cl_build_options.str().c_str()); 
    std::cout << "building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << "\n";
    if (cl_err!= CL_SUCCESS)
    {
        std::cout << " Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << "\n";
        exit(1);
    }

    context_g = context;
    default_device_g=  default_device;
    program_g= program;
}
