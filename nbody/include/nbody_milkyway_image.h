
#ifndef _NBODY_MILKYWAY_IMAGE_H_
#define _NBODY_MILKYWAY_IMAGE_H_

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

extern const unsigned char* milkywayImageRLEPixelData;

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

