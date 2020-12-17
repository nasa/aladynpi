c------------------------------------------------------------------
c 11-20-2020
c
c System Module Unit for aladyn_pi.f code
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
c Use to overwrite some OpenMP functions when compiled for OpenACC
c------------------------------------------------------------------
c
      MODULE sys_OMP

      use omp_lib

      save

      CONTAINS
c
c------------------------------------------------------------------
c
      integer function GET_NUM_PROCS()

       GET_NUM_PROCS = omp_get_num_procs()
       return

      end function 
c
c------------------------------------------------------------------
c
      integer function GET_NUM_THREADS()

       GET_NUM_THREADS = omp_get_num_threads()
       return

      end function 
c
c------------------------------------------------------------------
c
      integer function GET_MAX_THREADS()

       GET_MAX_THREADS = omp_get_max_threads()
       return

      end function 
c
c------------------------------------------------------------------
c
      subroutine SET_NUM_THREADS(num_thrds)

       integer, intent(in) :: num_thrds

       call omp_set_num_threads(num_thrds)
       return

      end subroutine
c
c------------------------------------------------------------------
c
      integer function GET_THREAD_NUM()

       GET_THREAD_NUM = OMP_get_thread_num()
       return

      end function 
c
c------------------------------------------------------------------
c
      double precision function GET_WTIME()

       GET_WTIME = OMP_get_wtime()
       return

      end function 
c------------------------------------------------------------------
c
      END MODULE  ! sys_OMP !
c
c =====================================================================
c
      MODULE sys_ACC

      use omp_lib
      use openacc
      use iso_c_binding

      save

      integer(kind=acc_device_kind) :: devicetype
      integer(kind=8) :: My_GPU_mem,My_GPU_free_mem

      CONTAINS
c
c------------------------------------------------------------------
c
      integer function gethostid() BIND(C)
       use iso_c_binding
       integer (C_INT) :: gethostid
      end function gethostid
c
c------------------------------------------------------------------
c
      integer function get_device_type()

       get_device_type = acc_get_device_type()
       return

      end function 
c
c------------------------------------------------------------------
c
      integer function get_num_devices(idev_type)

       get_num_devices = acc_get_num_devices(devicetype)

       return
      end function 
c
c------------------------------------------------------------------
c
      integer(kind=8) function get_gpu_mem(My_GPU)

       integer, intent(in) :: My_GPU

       get_gpu_mem = ACC_GET_PROPERTY(My_GPU, acc_device_current,
     1               acc_property_memory)
       return
      end function 
c
c------------------------------------------------------------------
c
      integer(kind=8)function get_gpu_free_mem(My_GPU)

       integer, intent(in) :: My_GPU

       get_gpu_free_mem =  ACC_GET_PROPERTY(My_GPU, acc_device_current,
     1                     acc_property_free_memory)
       return
      end function 
c
c------------------------------------------------------------------
c
      subroutine set_device_num(my_device, dev_type)

      integer, intent(in) :: my_device
      integer(kind=acc_device_kind), intent(in) :: dev_type

       call ACC_SET_DEVICE_NUM(my_device, dev_type)

      return
      end subroutine
c
c------------------------------------------------------------------
c
      subroutine GPU_Init(my_device)

       call ACC_INIT(my_device)

      return
      end subroutine
c
c------------------------------------------------------------------
c
      END MODULE  ! sys_ACC !
c
c =====================================================================
c
