
wedge = 82

start = {
   background = {
      q = 0.6,
      r0 = 20
   },
   streams = {
      {
         epsilon = -1,
         mu      = 31,
         r       = 29,
         theta   = 2,
         phi     = 0.1,
         sigma   = 2.8
      }
   }
}

max = {
   background = {
      q = 1.0000000000,
      r0 = 30.0000000000
   },
   streams = {
      {
         epsilon = 20.000000,
         mu      = 76.0000000000,
         r       = 45.7000000000,
         theta   = 6.2831853072,
         phi     = 6.2831853072,
         sigma   = 20.0000000000
      }
   }
}

min = {
   background = {
      q = 0.3000000000,
      r0 = 1.0000000000
   },
   streams = {
      {
         epsilon = -20.000000,
         mu      = -53.0000000000,
         r       = 0.7000000000,
         theta   = 0.0000000000,
         phi     = 0.0000000000,
         sigma   = 1.0000000000
      }
   }
}

step = {
   background = {
      q = 0.1000000,
      r0 = 0.1000000
   },
   streams = {
      {
         epsilon = 0.1000,
         mu      = 0.30000000,
         r       = 0.30000000,
         theta   = 0.30000000,
         phi     = 0.30000000,
         sigma   = 0.30000000
      }
   }
}

area = {
   {
      r_min = 16.0,
      r_max = 22.5,
      r_steps = 140,

      mu_min = 310,
      mu_max = 419,
      mu_steps = 160,

      nu_min = -1.25,
      nu_max = 1.25,
      nu_steps = 64
   }
}