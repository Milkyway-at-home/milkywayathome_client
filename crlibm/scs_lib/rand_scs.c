/*
 * Author  : Defour David
 * Contact : David.Defour@ens-lyon.fr
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
 */
#include <stdlib.h>
#include "scs.h"
#include "scs_private.h"

/*
 * Return 'sizeof(int)' random bits   
 */


int rand_val(void){
  int val;
  int i;

  val = (rand() & 0x000000ff);
  for(i=0; i<(sizeof(int)); i++){
    val = val << 8;
    val += (rand() & 0x000000ff ); /* we keep only 8 bits */
  }
  return val;
}


/*
 * Put into 'result' a scs random number with the index field set
 * with a value between -expo_max and +expo_max.
 *
 * Rem. :
 * 1) If you want an scs number belonging to double precision floating
 * point number you must call scs_rand with an expo_max less than 39.
 * 2) expo_max must be less than RAND_MAX that is usually set a
 * value greater than 32767
 */
void scs_rand(scs_ptr result, int expo_max){
  int i;

  R_EXP = 1;
  R_IND = (rand() % (2*expo_max)) - expo_max;
  R_SGN = ((2*rand()- RAND_MAX) > 0) ?  (-1) : (1);

  
  for(i=0; i<SCS_NB_WORDS; i++){
    /* We keep the first SCS_NB_BITS bits of a random value */
    R_HW[i] = rand_val() & SCS_MASK_RADIX; 
  }
}
