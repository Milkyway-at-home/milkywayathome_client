#include <CL/cl.h>
#include <CL/cl_ext.h>

#include "../astronomy/parameters.h"
#include "../astronomy/star_points.h"
#include "../astronomy_gpu/pi_constants.h"
#include "../astronomy_gpu/r_constants.h"

#define STREAM_PARAMETERS(x, y) (parameters[((x) * 6) + (y) + 3])
#define STREAM_WEIGHTS(x) parameters[((x) * 6) + 2]
#define BACKGROUND_WEIGHT (0.0)

#define LOCAL_WORK_SIZE 64

#define NVIDIA_INTEGRAL_KERNEL "nvidia_integral_kernel.cl"
#define NVIDIA_INTEGRAL_KERNEL2 "nvidia_integral_kernel2.cl"
#define ATI_INTEGRAL_KERNEL "ati_integral_kernel2.cl"
#define INTEGRAL_KERNEL_NAME "integral_kernel"

#define NVIDIA_LIKELIHOOD_KERNEL "nvidia_likelihood_kernel.cl"
#define ATI_LIKELIHOOD_KERNEL "ati_likelihood_kernel.cl"
#define LIKELIHOOD_KERNEL_NAME "likelihood_kernel"

#define ZERO_INTEGRAL_KERNEL "zero_integral_kernel.cl"
#define ZERO_INTEGRAL_KERNEL_NAME "zero_integral_kernel"

#define ZERO_LIKELIHOOD_KERNEL "zero_likelihood_kernel.cl"
#define ZERO_LIKELIHOOD_KERNEL_NAME "zero_likelihood_kernel"

enum PLATFORM {
  NVIDIA,
  ATI
};

typedef struct {
  //integral memory
  cl_mem *sinb, *sinl, *cosb, *cosl;
  cl_mem *v;
  cl_mem *r_point, *qw_r3_N;
  cl_mem *bg_int, *st_int;
  cl_mem fstream_a, fstream_c, inv_fstream_sigma_sq2;
  //likelihood memory
  cl_mem stars;
  cl_mem qgaus_W, dx;
  cl_mem probability;
  cl_mem bg_weight, st_weight;
  cl_mem A;
  //book keeping
  cl_context context;
  cl_command_queue queue;
  cl_device_id *devices;
  int num_integrals;
  PLATFORM platform;
} ocl_mem_t;

cl_platform_id * get_platforms();

char * get_platform_info(cl_platform_id platform,
			 cl_platform_info info);

cl_context get_context(cl_device_id device);

cl_device_id * get_devices(cl_platform_id platform);

void setup_command_queue(ocl_mem_t *ocl_mem);

ocl_mem_t * setup_ocl(ASTRONOMY_PARAMETERS *ap,
		      STAR_POINTS *sp);

void setup_ocl(double *parameters, ocl_mem_t *ocl_mem,
	       int number_streams, int sgr_coordinates,
	       int wedge);

void setup_weights(ocl_mem_t *ocl_mem,
		   double *parameters,
		   double *bg_int,
		   double *st_int,
		   int number_streams);

void destruct_ocl(ocl_mem_t *ocl_mem);

double ocl_likelihood(double *parameters,
		      ASTRONOMY_PARAMETERS *ap,
		      STAR_POINTS *sp,
		      ocl_mem_t *ocl_mem);

const char * read_kernel(const char *kernel_source);

void build_kernel(cl_program program, 
		  cl_device_id *devices,
		  const char *name);

int determine_work_size(int local_work_size, int desired_size);

#define check_error(err)					\
{								\
  if (err != CL_SUCCESS)					\
    {								\
      fprintf(stdout, "Error in OpenCL %d at %s:%s:%d\n", (err),	\
	      __FILE__, __FUNCTION__, __LINE__);			\
      exit(1);								\
    }									\
}

