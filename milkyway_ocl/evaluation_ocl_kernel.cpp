/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string>
#include <fstream>
#include "boinc_api.h"

#include "evaluation_ocl.h"

const char * read_kernel(const char *kernel_source)
{
#ifdef BOINC_APPLICATION
  std::string real_path;
  int ret = boinc_resolve_filename_s(kernel_source, real_path);
  if (ret)
    {
      fprintf(stderr, "Unable to resolve the path of %s\n",
	      kernel_source);
      return 0;
    }
  else
#else
    std::string real_path = kernel_source;
#endif
    {
      FILE *fp = fopen(real_path.c_str(), "rb");
      if (fp)
	{
	  fseek(fp, 0, SEEK_END);
	  unsigned int size = ftell(fp);
	  rewind(fp);
	  char *buffer = new char[size+1];
	  fread((void*)buffer, size, size, fp);
	  fclose(fp);
	  buffer[size] = '\0';
	  return buffer;
	}
      else
	{
	  fprintf(stderr, "Error opening OpenCL kernel %s\n",
		  real_path.c_str());
	  return 0;
	}
    }
}

void write_binaries(cl_program program,
		    const char *name)
{
  cl_uint num_devices;
  clGetProgramInfo(program, CL_PROGRAM_NUM_DEVICES,
		   sizeof(cl_uint), &num_devices, 0);
  size_t *binary_size = new size_t[num_devices];
  clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES,
		   num_devices * sizeof(size_t), binary_size, 0);
  unsigned char **binaries = new unsigned char*[num_devices];
  size_t total_size;
  for(unsigned int  i = 0;i<num_devices;++i)
    {
      printf("Binary %d is %u bytes\n", i, (unsigned int) binary_size[i]);
      binaries[i] = new unsigned char[binary_size[i]];
      total_size += binary_size[i];
    }
  clGetProgramInfo(program, CL_PROGRAM_BINARIES,
		   total_size, binaries, 0);
  for(unsigned int  i = 0;i<num_devices;++i)
    {
      std::ofstream output(name);
      if (output)
	{
	  output.write((const char*) binaries[i], binary_size[i]);
	}
      delete [] binaries[i];
    }
  delete [] binary_size;
  delete [] binaries;
}

void build_kernel(cl_program program,
		  cl_device_id *devices,
		  const char *name)
{
  // char build_options[] = "-cl-single-precision-constant "
  //   "-cl-strict-aliasing "
  //   "-cl-mad-enable "
  //   "-cl-unsafe-math-optimizations "
  //   "-cl-finite-math-only "
  //   "-cl-fast-relaxed-math "
  //   "-Werror ";
  char build_options[] = "-w -I .";
  cl_int err = clBuildProgram(program, 1,
			      devices, build_options, 0, 0);
  //check the build log
  if (err != CL_SUCCESS)
    printf("Error Building Kernel, error code %d\n", err);
  size_t len;
  clGetProgramBuildInfo(program,
			devices[0],
			CL_PROGRAM_BUILD_LOG,
			0, 0, &len);
  char *error_message = new char[len];
  clGetProgramBuildInfo(program,
			devices[0],
			CL_PROGRAM_BUILD_LOG,
			len, error_message, 0);
  printf("%s\n", error_message);
  delete [] error_message;
  if (err != CL_SUCCESS)
    exit(1);
  //write_binaries(program, name);
}

