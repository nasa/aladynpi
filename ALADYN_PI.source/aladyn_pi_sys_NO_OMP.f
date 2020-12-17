c------------------------------------------------------------------
c 11-20-2020
c
c System Module Unit for the aladyn_pi.f code
c
c Vesselin Yamakov
c National Institute of Aerospace
c 100 Exploration Way,
c Hampton, VA 23666 
c phone: (757)-864-2850
c fax:   (757)-864-8911
c e-mail: yamakov@nianet.org
c
c------------------------------------------------------------------
c Use to overwrite some OpenMP and OpenACC functions when not available
c------------------------------------------------------------------
c
      MODULE sys_OMP

      save

      CONTAINS
c
c------------------------------------------------------------------
c
      integer function GET_NUM_PROCS()

       GET_NUM_PROCS = 1
       return

      end function 
c
c------------------------------------------------------------------
c
      integer function GET_NUM_THREADS()

       GET_NUM_THREADS = 1
       return

      end function 
c
c------------------------------------------------------------------
c
      integer function GET_MAX_THREADS()

       GET_MAX_THREADS = 1
       return

      end function 
c
c------------------------------------------------------------------
c
      subroutine SET_NUM_THREADS(num_thrds)

       integer, intent(in) :: num_thrds

       return

      end subroutine
c
c------------------------------------------------------------------
c
      integer function GET_THREAD_NUM()

       GET_THREAD_NUM = 0
       return

      end function 
c
c------------------------------------------------------------------
c
      double precision function GET_WTIME()

       call cpu_time(GET_WTIME)
       return

      end function 
c
c------------------------------------------------------------------
c
      END MODULE  ! sys_OMP !
c
c =====================================================================
c Overwrites some OpenACC functions when not available
c------------------------------------------------------------------
c

      MODULE sys_ACC

      save

      integer :: devicetype = 0
      integer :: acc_device_current = 0
      integer(kind=8) :: My_GPU_mem,My_GPU_free_mem

      CONTAINS
c
c------------------------------------------------------------------
c
      integer function get_device_type()

       get_device_type = 0
       return

      end function 
c
c------------------------------------------------------------------
c
      integer function get_num_devices(idev_type)

       get_num_devices = 0
       return

      end function 
c
c------------------------------------------------------------------
c
      subroutine set_device_num(my_device, dev_type)

      integer, intent(in) :: my_device, dev_type

      return
      end subroutine
c
c------------------------------------------------------------------
c
      integer(kind=8) function get_gpu_mem(My_GPU)
       integer, intent(in) :: My_GPU

        get_gpu_mem = 0

       return
      end function 
c
c------------------------------------------------------------------
c
      integer(kind=8)function get_gpu_free_mem(My_GPU)
       integer, intent(in) :: My_GPU

        get_gpu_free_mem = 0

       return
      end function 
c
c------------------------------------------------------------------
c
      subroutine GPU_Init(my_device)
       integer, intent(in) :: my_device

       My_GPU_mem = 0; My_GPU_free_mem = 0

      return
      end subroutine
c
c------------------------------------------------------------------
c
      END MODULE  ! sys_ACC !
c
c =====================================================================
c
