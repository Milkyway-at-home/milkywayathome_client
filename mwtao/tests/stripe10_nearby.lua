
wedge = 10

start = {
   background = {
      q  = 0.5,
      r0 = 8.0
   },
   streams = {
      {
         epsilon = -1.3,
         mu      = 230.0,
         r       = 36.0,
         theta   = -1.5,
         phi     = -0.1,
         sigma   = 5.0
      },

      {
         epsilon = -2.0,
         mu      = 195.0,
         r       = 10.0,
         theta   = 0.15,
         phi     = -3.5,
         sigma   = 22.0
      },

      {
         epsilon = -2.0,
         mu      = 186.0,
         r       = 20,
         theta   = 0.3,
         phi     = 2.0,
         sigma   = 6.0
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
         mu      = 227.0,
         r       = 72.4,
         theta   = 6.283185307179586,
         phi     = 6.283185307179586,
         sigma   = 25.0
      },

      {
         epsilon = 20.0,
         mu      = 227.0,
         r       = 72.4,
         theta   = 6.283185307179586,
         phi     = 6.283185307179586,
         sigma   = 25.0
      },

      {
         epsilon = 20.0,
         mu      = 227.0,
         r       = 72.4,
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
         mu      = 165.0,
         r       = 2.30,
         theta   = -6.283185307179586,
         phi     = -6.283185307179586,
         sigma   = 0.10
      },

      {
         epsilon = -20.0,
         mu      = 165.0,
         r       = 2.30,
         theta   = -6.283185307179586,
         phi     = -6.283185307179586,
         sigma   = 0.10
      },

      {
         epsilon = -20.0,
         mu      = 165.0,
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
      r_max = 23.0,
      r_steps = 175,

      mu_min = 165,
      mu_max = 245,
      mu_steps = 400,

      nu_min = -1.25,
      nu_max = 1.25,
      nu_steps = 40
   },

   {
      r_min = 16.0,
      r_max = 23.0,
      r_steps = 95,

      mu_min = 227.0,
      mu_max = 230.0,
      mu_steps = 200,

      nu_min = -1.25,
      nu_max = 0.5,
      nu_steps = 20
   }
}
