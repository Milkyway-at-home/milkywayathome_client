
wedge = 22

start = {
   background = {
      q  = 0.3,
      r0 = 2.9
   },
   streams = {
      {
         epsilon = -1.7,
         mu      = 136.0,
         r       = 18.0,
         theta   = 0.72,
         phi     = -1.24,
         sigma   = 1.7
      },

      {
         epsilon = 0.6,
         mu      = 184.0,
         r       = -3.0,
         theta   = 1.15,
         phi     = 2.1,
         sigma   = 24.0
      },

      {
         epsilon = 0.8,
         mu      = 132.0,
         r       = -26.0,
         theta   = 0.85,
         phi     = -2.2,
         sigma   = 7.0
      }
   }
}

max = {
   background = {
      q  = 1.0,
      r0 = 30.0
   },
   streams = {
      {
         epsilon = 20.0,
         mu      = 249.0,
         r       = 45.6,
         theta   = 6.283185307179586,
         phi     = 6.283185307179586,
         sigma   = 25.0
      },

      {
         epsilon = 20.0,
         mu      = 249.0,
         r       = 45.6,
         theta   = 6.283185307179586,
         phi     = 6.283185307179586,
         sigma   = 25.0
      },

      {
         epsilon = 20.0,
         mu      = 249.0,
         r       = 45.6,
         theta   = 6.283185307179586,
         phi     = 6.283185307179586,
         sigma   = 25.0
      }
   }
}

min = {
   background = {
      q  = 0.3,
      r0 = 1.0
   },
   streams = {
      {
         epsilon = -20.0,
         mu      = 131.0,
         r       = 2.30,
         theta   = -6.283185307179586,
         phi     = -6.283185307179586,
         sigma   = 0.10
      },

      {
         epsilon = -20.0,
         mu      = 131.0,
         r       = 2.30,
         theta   = -6.283185307179586,
         phi     = -6.283185307179586,
         sigma   = 0.10
      },

      {
         epsilon = -20.0,
         mu      = 131.0,
         r       = 2.30,
         theta   = -6.283185307179586,
         phi     = -6.283185307179586,
         sigma   = 0.10
      }
   }
}


step = {
   background = {
      q  = 0.000004,
      r0 = 0.00008
   },
   streams = {
      {
         epsilon = 0.000001,
         mu      = 0.00003,
         r       = 0.00004,
         theta   = 0.00006,
         phi     = 0.00004,
         sigma   = 0.000004
      },

      {
         epsilon = 0.000001,
         mu      = 0.00003,
         r       = 0.00004,
         theta   = 0.00006,
         phi     = 0.00004,
         sigma   = 0.000004
      },

      {
         epsilon = 0.000001,
         mu      = 0.00003,
         r       = 0.00004,
         theta   = 0.00006,
         phi     = 0.00004,
         sigma   = 0.000004
      }
   }
}

area = {
   {
      r_min = 16.0,
      r_max = 22.5,
      r_steps = 140,

      mu_min = 131,
      mu_max = 225,
      mu_steps = 160,

      nu_min = -1.25,
      nu_max = 1.25,
      nu_steps = 64
   },

   {
      r_min = 16.0,
      r_max = 22.5,
      r_steps = 140,

      mu_min = 207.0,
      mu_max = 209.0,
      mu_steps = 40,

      nu_min = 0.8,
      nu_max = 1.25,
      nu_steps = 32
   },

   {
      r_min = 16.0,
      r_max = 22.5,
      r_steps = 140,

      mu_min = 202.0,
      mu_max = 204.0,
      mu_steps = 40,

      nu_min = -0.5,
      nu_max = 0.4,
      nu_steps = 32
   }
}
