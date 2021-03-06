*
* samplexy_curv.F
*
* Ansley Manke
* 4/03 
* ???  Note: does not currently deal with modulo data  ???
* 1/06 Remove calls to SIGN to compile with f90 on IRIX.

* This function samples 4-d curvlinear data on the x and y axes 
* indicated by args 4 and 5
*
* In this subroutine we provide information about
* the function.  The user configurable information 
* consists of the following:
*
* descr              Text description of the function
*
* num_args           Required number of arguments
*
* axis_inheritance   Type of axis for the result
*                       ( CUSTOM, IMPLIED_BY_ARGS, NORMAL, ABSTRACT )
*                       CUSTOM          - user defined axis
*                       IMPLIED_BY_ARGS - same axis as the incoming argument
*                       NORMAL          - the result is normal to this axis
*                       ABSTRACT        - an axis which only has index values
*
* piecemeal_ok       For memory optimization:
*                       axes where calculation may be performed piecemeal
*                       ( YES, NO )
* 
*
* For each argument we provide the following information:
*
* name               Text name for an argument
*
* unit               Text units for an argument
*
* desc               Text description of an argument
*
* axis_influence     Are this argument's axes the same as the result grid?
*                       ( YES, NO )
*
* axis_extend       How much does Ferret need to extend arg limits relative to result 
*


      SUBROUTINE samplexy_curv_init(id)

      INCLUDE 'ferret_cmn/EF_Util.cmn'

      INTEGER id, arg

***********************************************************************
*                                           USER CONFIGURABLE PORTION |
*                                                                     |
*                                                                     V
      CHARACTER*100 fcn_desc
      WRITE (fcn_desc, 10)
   10 FORMAT ('Returns data sampled at a set of (X,Y) points, ',
     .   'using linear interpolation')
      CALL ef_set_desc(id, fcn_desc)

      CALL ef_set_num_args(id, 5)
      CALL ef_set_has_vari_args(id, NO)
      CALL ef_set_axis_inheritance(id, ABSTRACT, 
     .     NORMAL, IMPLIED_BY_ARGS, IMPLIED_BY_ARGS)

      CALL ef_set_num_work_arrays(id, 1)

      CALL ef_set_piecemeal_ok(id, NO, NO, NO, NO)

      arg = 1
      CALL ef_set_arg_name(id, arg, 'DAT_TO_SAMPLE')
      CALL ef_set_arg_desc(id, arg, 'variable (x,y,z,t) to sample')
      CALL ef_set_axis_influence(id, arg, NO, NO, YES, YES)

      arg = 2
      CALL ef_set_arg_name(id, arg, 'DAT_LON')
      CALL ef_set_arg_desc(id, arg, 'longitudes of input variable')
      CALL ef_set_axis_influence(id, arg, NO, NO, YES, YES)

      arg = 3
      CALL ef_set_arg_name(id, arg, 'DAT_LAT')
      CALL ef_set_arg_desc(id, arg, 'latitudes of input variable')
      CALL ef_set_axis_influence(id, arg, NO, NO, YES, YES)

      arg = 4
      CALL ef_set_arg_name(id, arg, 'XPTS')
      CALL ef_set_arg_desc(id, arg, 'X values of sample points')
      CALL ef_set_axis_influence(id, arg, NO, NO, NO, NO)

      arg = 5
      CALL ef_set_arg_name(id, arg, 'YPTS')
      CALL ef_set_arg_desc(id, arg, 'Y values of sample points')
      CALL ef_set_axis_influence(id, arg, NO, NO, NO, NO)
*                                                                     ^
*                                                                     |
*                                           USER CONFIGURABLE PORTION |
***********************************************************************

      RETURN 
      END

*
* In this subroutine we provide information about the lo and hi
* limits associated with each abstract or custom axis.   The user 
* configurable information consists of the following:
*
* loss               lo subscript for an axis
*
* hiss               hi subscript for an axis
*

      SUBROUTINE samplexy_curv_result_limits(id)

      INCLUDE 'ferret_cmn/EF_Util.cmn'

      INTEGER id
      INTEGER arg_lo_ss(4,EF_MAX_ARGS), arg_hi_ss(4,EF_MAX_ARGS),
     .     arg_incr(4,EF_MAX_ARGS)

* **********************************************************************
*                                           USER CONFIGURABLE PORTION |
*                                                                     |
*                                                                     V

      INTEGER my_lo_l, my_hi_l
      INTEGER nx, ny, nz, nt

*     Use utility functions to get context information about the 
*     1st argument, to set the abstract axis lo and hi indices.

      CALL ef_get_arg_subscripts(id, arg_lo_ss, arg_hi_ss, arg_incr)

      nx = arg_hi_ss(X_AXIS, ARG4) - arg_lo_ss(X_AXIS, ARG4) + 1
      ny = arg_hi_ss(Y_AXIS, ARG4) - arg_lo_ss(Y_AXIS, ARG4) + 1
      nz = arg_hi_ss(Z_AXIS, ARG4) - arg_lo_ss(Z_AXIS, ARG4) + 1
      nt = arg_hi_ss(T_AXIS, ARG4) - arg_lo_ss(T_AXIS, ARG4) + 1

      my_lo_l = 1
      my_hi_l = max(nx,ny,nz,nt)

      CALL ef_set_axis_limits(id, X_AXIS, my_lo_l, my_hi_l)

*                                                                     ^
*                                                                     |
*                                           USER CONFIGURABLE PORTION |
* **********************************************************************

      RETURN 
      END

*
* In this subroutine we request an amount of storage to be supplied
* by Ferret and passed as an additional argument.
*
      SUBROUTINE samplexy_curv_work_size(id)

      INCLUDE 'ferret_cmn/EF_Util.cmn'
      INCLUDE 'ferret_cmn/EF_mem_subsc.cmn'

      INTEGER id

* **********************************************************************
*                                            USER CONFIGURABLE PORTION |
*                                                                      |
*                                                                      V

*
* Set the work arrays,  X/Y/Z/T dimensions
*
* ef_set_work_array_lens(id,array #,xlo,ylo,zlo,tlo,xhi,yhi,zhi,thi)
*
      INTEGER nx, ny

      INTEGER arg_lo_ss(4,1:EF_MAX_ARGS), arg_hi_ss(4,1:EF_MAX_ARGS),
     .     arg_incr(4,1:EF_MAX_ARGS)

      CALL ef_get_arg_subscripts(id, arg_lo_ss, arg_hi_ss, arg_incr)

*  workspace for saving locations of xpt in arg2 array.

      nx = arg_hi_ss(X_AXIS,ARG2) - arg_lo_ss(X_AXIS,ARG2) 
      ny = arg_hi_ss(Y_AXIS,ARG2) - arg_lo_ss(Y_AXIS,ARG2) 

*  xlocations
      CALL ef_set_work_array_dims (id, 1, 1, 1, 1, 1, nx*ny, 2, 1, 1)

*                                                                      ^
*                                                                      |
*                                            USER CONFIGURABLE PORTION |
* **********************************************************************

      RETURN 
      END
*
* In this subroutine we compute the result
*



      SUBROUTINE samplexy_curv_compute(id, arg_1, arg_2, arg_3, arg_4, 
     .       arg_5, result, xlocations)

      INCLUDE 'ferret_cmn/EF_Util.cmn'
      INCLUDE 'ferret_cmn/EF_mem_subsc.cmn'

      INTEGER id

      REAL bad_flag(EF_MAX_ARGS), bad_flag_result
      REAL arg_1(mem1lox:mem1hix, mem1loy:mem1hiy, mem1loz:mem1hiz,
     .     mem1lot:mem1hit)
      REAL arg_2(mem2lox:mem2hix, mem2loy:mem2hiy, mem2loz:mem2hiz,
     .     mem2lot:mem2hit)
      REAL arg_3(mem3lox:mem3hix, mem3loy:mem3hiy, mem3loz:mem3hiz,
     .     mem3lot:mem3hit)
      REAL arg_4(mem4lox:mem4hix, mem4loy:mem4hiy, mem4loz:mem4hiz,
     .     mem4lot:mem4hit)
      REAL arg_5(mem5lox:mem5hix, mem5loy:mem5hiy, mem5loz:mem5hiz,
     .     mem5lot:mem5hit)
      REAL result(memreslox:memreshix, memresloy:memreshiy,
     .      memresloz:memreshiz, memreslot:memreshit)


* After initialization, the 'res_' arrays contain indexing information
* for the result axes.  The 'arg_' arrays will contain the indexing
* information for each variable's axes.

      INTEGER res_lo_ss(4), res_hi_ss(4), res_incr(4)
      INTEGER arg_lo_ss(4,EF_MAX_ARGS), arg_hi_ss(4,EF_MAX_ARGS),
     .     arg_incr(4,EF_MAX_ARGS)

***********************************************************************
*                                           USER CONFIGURABLE PORTION |
*                                                                     |
*                                                                     V
      COMMON /STOR/ mxdat, mydat
      INTEGER mxdat, mydat
      INTEGER nx, nxx, nxy, nxz, nxt
      INTEGER ny, nyx, nyy, nyz, nyt
      INTEGER ndimx, ndimy
      INTEGER nx1, nx2, nx3, ny1, ny2, ny3

*  Set up work arrays

      REAL xlocations(wrk1lox:wrk1lox+(wrk1hix-wrk1lox),wrk1loy:wrk1hiy,
     .                wrk1loz:wrk1hiz, wrk1lot:wrk1hit)
      REAL ylocations(2)

      INTEGER i, j, k, l
      INTEGER i1,j1,k1,l1
      INTEGER i2,j2,k2,l2
      INTEGER i3,j3,k3,l3
      INTEGER i4,j4,k4,l4
      INTEGER i5,j5,k5,l5

      INTEGER ibot, itop, jbot, jtop
      INTEGER nloc, iloc, ptloc, icross
      REAL xpt, ypt, ydiff, yd
      REAL fxbot, fxtop, fbb, ftb, fbt, ftt
      REAL xbot, xtop, ybot, ytop
      REAL frac
      CHARACTER*255 err_msg

      CALL ef_get_res_subscripts(id, res_lo_ss, res_hi_ss, res_incr)
      CALL ef_get_arg_subscripts(id, arg_lo_ss, arg_hi_ss, arg_incr)
      CALL ef_get_bad_flags(id, bad_flag, bad_flag_result)

C Check that first 3 args are on the same grid.

      nx1 = arg_hi_ss(X_AXIS,ARG1) - arg_lo_ss(X_AXIS,ARG1) + 1
      nx2 = arg_hi_ss(X_AXIS,ARG2) - arg_lo_ss(X_AXIS,ARG2) + 1
      nx3 = arg_hi_ss(X_AXIS,ARG3) - arg_lo_ss(X_AXIS,ARG3) + 1
      ny1 = arg_hi_ss(X_AXIS,ARG1) - arg_lo_ss(X_AXIS,ARG1) + 1
      ny2 = arg_hi_ss(X_AXIS,ARG2) - arg_lo_ss(X_AXIS,ARG2) + 1
      ny3 = arg_hi_ss(X_AXIS,ARG3) - arg_lo_ss(X_AXIS,ARG3) + 1

      IF (nx1 .NE. nx2 .OR. nx2 .NE. nx3 .OR. nx1 .NE. nx3) THEN
         err_msg ='Arguments 1, 2, and 3 must have the same XY grid'
         GO TO 999
      ENDIF
      IF (ny1 .NE. ny2 .OR. ny2 .NE. ny3 .OR. ny1 .NE. ny3) THEN
         err_msg ='Arguments 1, 2, and 3 must have the same XY grid'
         GO TO 999
      ENDIF

C Check that sample x and y are simple lists of points - same length

      nxx = arg_hi_ss(X_AXIS,ARG4) - arg_lo_ss(X_AXIS,ARG4) + 1
      nxy = arg_hi_ss(Y_AXIS,ARG4) - arg_lo_ss(Y_AXIS,ARG4) + 1
      nxz = arg_hi_ss(Z_AXIS,ARG4) - arg_lo_ss(Z_AXIS,ARG4) + 1
      nxt = arg_hi_ss(T_AXIS,ARG4) - arg_lo_ss(T_AXIS,ARG4) + 1
      nx = max(nxx, nxy, nxz, nxt)

      nyx = arg_hi_ss(X_AXIS,ARG5) - arg_lo_ss(X_AXIS,ARG5) + 1
      nyy = arg_hi_ss(Y_AXIS,ARG5) - arg_lo_ss(Y_AXIS,ARG5) + 1
      nyz = arg_hi_ss(Z_AXIS,ARG5) - arg_lo_ss(Z_AXIS,ARG5) + 1
      nyt = arg_hi_ss(T_AXIS,ARG5) - arg_lo_ss(T_AXIS,ARG5) + 1
      ny = max(nyx, nyy, nyz, nyt)

      ndimx = 0
      ndimy = 0
      DO 110 i = X_AXIS,T_AXIS
         IF(arg_hi_ss(i,ARG4) - arg_lo_ss(i,ARG4) .GT.0) ndimx=ndimx + 1
         IF(arg_hi_ss(i,ARG5) - arg_lo_ss(i,ARG5) .GT.0) ndimy=ndimy + 1
  110 CONTINUE

      IF (nx .NE. ny   .OR.  ndimx .GT. 1  .OR. ndimy .GT.1) THEN
         err_msg = 
     .  'Arguments 4 and 5 must be 1-dimensional lists of equal length'
         GO TO 999
      ENDIF


*  For each (xpt,ypt) pair, search the data arrays arg2,arg3
*  for the nearest higher (x,y) lat-lon coordinates.  Interpolate 
*  arg_1 in 2 directions for the result.

      i4 = arg_lo_ss(X_AXIS,ARG4)
      j4 = arg_lo_ss(Y_AXIS,ARG4)
      k4 = arg_lo_ss(Z_AXIS,ARG4)
      l4 = arg_lo_ss(T_AXIS,ARG4)

      i5 = arg_lo_ss(X_AXIS,ARG5)
      j5 = arg_lo_ss(Y_AXIS,ARG5)
      k5 = arg_lo_ss(Z_AXIS,ARG5)
      l5 = arg_lo_ss(T_AXIS,ARG5)

      j = res_lo_ss(Y_AXIS)
      DO 600 i = res_lo_ss(X_AXIS), res_hi_ss(X_AXIS)

* Find all locations of xpt in arg_2, longitudes of data

	 xpt = arg_4(i4,j4,k4,l4)

         k2 = arg_lo_ss(Z_AXIS,ARG2)
         l2 = arg_lo_ss(T_AXIS,ARG2)

         nloc = 0

         DO 200 j2 = arg_lo_ss(Y_AXIS,ARG2), arg_hi_ss(Y_AXIS,ARG2)
            i2 = arg_lo_ss(X_AXIS,ARG2)
c            icross = SIGN(1, (xpt - arg_2(i2,j2,k2,l2)) )
            icross = 1
            IF ((xpt - arg_2(i2,j2,k2,l2)) .LT. 0.) icross = -1

            DO 100 i2 = arg_lo_ss(X_AXIS,ARG2), arg_hi_ss(X_AXIS,ARG2)
c               iloc = SIGN(1, (xpt-arg_2(i2,j2,k2,l2)) )
               iloc = 1
               IF ((xpt-arg_2(i2,j2,k2,l2)) .LT. 0. ) iloc = -1

               IF (icross .NE. iloc)THEN
                  nloc = nloc + 1
                  xlocations(nloc,1,1,1) = MAX(i2-1,1)
                  xlocations(nloc,2,1,1) = j2
c                  icross = SIGN(1, (xpt - arg_2(i2,j2,k2,l2)) )
                  icross = 1
                  IF ((xpt - arg_2(i2,j2,k2,l2)) .LT. 0.) icross = -1
               ENDIF

  100       CONTINUE
  200    CONTINUE

* Which of those locations is nearest xpt in arg_3, latitudes of data?

         ypt = arg_5(i5,j5,k5,l5)
         ydiff = ABS(bad_flag_result)
         ptloc = 1

         k3 = arg_lo_ss(Z_AXIS,ARG3)
         l3 = arg_lo_ss(T_AXIS,ARG3)

         DO 300 iloc = 1, nloc
            i3 = xlocations(iloc,1,1,1)
            j3 = xlocations(iloc,2,1,1)
               yd = ABS(ypt - arg_3(i3,j3,k3,l3))

            IF ( yd .LT. ydiff ) THEN
               ylocations(1) = i3
               ylocations(2) = j3
               ydiff = yd
               ptloc = iloc
            ENDIF
  300   CONTINUE

        ibot = xlocations(ptloc,1,1,1)
        jbot = xlocations(ptloc,2,1,1)
        IF (ibot .GE. arg_hi_ss(X_AXIS,ARG1) ) ibot = ibot - 1
        IF (jbot .GE. arg_hi_ss(Y_AXIS,ARG1) ) jbot = jbot - 1

        itop = ibot + 1
        jtop = jbot + 1

*  First interpolate in x, getting values of the fcn at (x,jbot)  and (x,jtop)


         k1 = arg_lo_ss(Z_AXIS,ARG1)
         k2 = arg_lo_ss(Z_AXIS,ARG2)
         k3 = arg_lo_ss(Z_AXIS,ARG3)
         DO 500 k = res_lo_ss(Z_AXIS), res_hi_ss(Z_AXIS)

            l1 = arg_lo_ss(T_AXIS,ARG1)
            l2 = arg_lo_ss(T_AXIS,ARG2)
            l3 = arg_lo_ss(T_AXIS,ARG3)
            DO 400 l = res_lo_ss(T_AXIS), res_hi_ss(T_AXIS)

               IF (ibot .GE. arg_lo_ss(X_AXIS,ARG1)  .AND. 
     .              itop .LT. arg_hi_ss(X_AXIS,ARG1) ) THEN

                  xbot = arg_2(ibot,jbot,k2,l2)
                  xtop = arg_2(itop,jbot,k2,l2)
      
                  fbb = arg_1(ibot,jbot,k1,l1) 
                  ftb = arg_1(itop,jbot,k1,l1)
                  fbt = arg_1(ibot,jtop,k1,l1)
                  ftt = arg_1(itop,jtop,k1,l1)

                  IF (fbb .NE. bad_flag(ARG1)  .AND.
     .                ftb .NE. bad_flag(ARG1)  .AND.
     .                fbt .NE. bad_flag(ARG1)  .AND.
     .                ftt .NE. bad_flag(ARG1) ) THEN


                     frac = (xpt - xbot ) / (xtop - xbot)
                     fxbot = fbb + frac* (ftb - fbb)
                     fxtop = fbt + frac* (ftt - fbt)

*  Now interpolate in y, getting value at (x,y)

                     IF (jbot .GE. arg_lo_ss(Y_AXIS,ARG1)  .AND. 
     .                   jtop .LE. arg_hi_ss(Y_AXIS,ARG1) ) THEN
                        ybot = arg_3(ibot,jbot,k3,l3)
                        ytop = arg_3(ibot,jtop,k3,l3)

                        frac = (ypt - ybot)/ (ytop - ybot)
                           result(i,j,k,l) = fxbot + frac* (fxtop-fxbot)
                        ELSE
                              result(i,j,k,l) = bad_flag_result
                        ENDIF

                     ELSE
                        result(i,j,k,l) = bad_flag_result
                     ENDIF		! bad_flag(ARG1) test
                  ENDIF


               l1 = l1 + arg_incr(T_AXIS,ARG1)
               l2 = l2 + arg_incr(T_AXIS,ARG2)
               l3 = l3 + arg_incr(T_AXIS,ARG3)

 400        CONTINUE

            k1 = k1 + arg_incr(Z_AXIS,ARG1)
            k2 = k2 + arg_incr(Z_AXIS,ARG2)
            k3 = k3 + arg_incr(Z_AXIS,ARG3)
 500     CONTINUE

         i4 = i4 + arg_incr(X_AXIS,ARG4)
         j4 = j4 + arg_incr(Y_AXIS,ARG4)
         k4 = k4 + arg_incr(Z_AXIS,ARG4)
         l4 = l4 + arg_incr(T_AXIS,ARG4)

         i5 = i5 + arg_incr(X_AXIS,ARG5)
         j5 = j5 + arg_incr(Y_AXIS,ARG5)
         k5 = k5 + arg_incr(Z_AXIS,ARG5)
         l5 = l5 + arg_incr(T_AXIS,ARG5)
600   CONTINUE

      RETURN
 999  CALL ef_bail_out (id, err_msg)

      END
