/*
  JPL Image Use Policy Unless otherwise noted, images and video on JPL
  public web sites (public sites ending with a jpl.nasa.gov address)
  may be used for any purpose without prior permission, subject to the
  special cases noted below. Publishers who wish to have authorization
  may print this page and retain it for their records; JPL does not
  issue image permissions on an image by image basis.

  By electing to download the material from this web site the user
  agrees:

  1. that Caltech makes no representations or warranties with respect
  to ownership of copyrights in the images, and does not represent
  others who may claim to be authors or owners of copyright of any of
  the images, and makes no warranties as to the quality of the
  images. Caltech shall not be responsible for any loss or expenses
  resulting from the use of the images, and you release and hold
  Caltech harmless from all liability arising from such use.

  2. to use a credit line in connection with images. Unless otherwise
  noted in the caption information for an image, the credit line
  should be "Courtesy NASA/JPL-Caltech."

  3. that the endorsement of any product or service by Caltech, JPL or
  NASA must not be claimed or implied.
*/

#ifndef _NBODY_MILKYWAY_IMAGE_H_
#define _NBODY_MILKYWAY_IMAGE_H_

/* GIMP RGB C-Source image dump 1-byte-run-length-encoded (milkyway_rle.c) */

#define MILKYWAYIMAGE_RUN_LENGTH_DECODE(image_buf, rle_data, size, bpp) do \
{ unsigned int __bpp; unsigned char *__ip; const unsigned char *__il, *__rd; \
  __bpp = (bpp); __ip = (image_buf); __il = __ip + (size) * __bpp; \
  __rd = (rle_data); if (__bpp > 3) { /* RGBA */ \
    while (__ip < __il) { unsigned int __l = *(__rd++); \
      if (__l & 128) { __l = __l - 128; \
        do { memcpy (__ip, __rd, 4); __ip += 4; } while (--__l); __rd += 4; \
      } else { __l *= 4; memcpy (__ip, __rd, __l); \
               __ip += __l; __rd += __l; } } \
  } else { /* RGB */ \
    while (__ip < __il) { unsigned int __l = *(__rd++); \
      if (__l & 128) { __l = __l - 128; \
        do { memcpy (__ip, __rd, 3); __ip += 3; } while (--__l); __rd += 3; \
      } else { __l *= 3; memcpy (__ip, __rd, __l); \
               __ip += __l; __rd += __l; } } \
  } } while (0)

#ifdef __cplusplus
extern "C" {
#endif

#define MILKYWAY_IMAGE_DATA_SIZE (2665007 + 1)

extern const unsigned char nbody_milkyway_image_bin[];

const unsigned char* milkywayImageRLEPixelData = nbody_milkyway_image_bin;

static const struct
{
  unsigned int 	 width;
  unsigned int 	 height;
  unsigned int 	 bytes_per_pixel; /* 3:RGB, 4:RGBA */
} milkywayImageParams = {
  1024, 1024, 3
};


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_MILKYWAY_IMAGE_H_ */

