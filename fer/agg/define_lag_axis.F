        SUBROUTINE DEFINE_LAG_AXIS(nsets, member_sets, memb_filenames,
     .                     fmatch, findices, tmatch, t1indices, tdims,
     .                     axcalf, axcalt, axlagf, axlagt, status)
*
*  This software was developed by the Thermal Modeling and Analysis
*  Project(TMAP) of the National Oceanographic and Atmospheric
*  Administration's (NOAA) Pacific Marine Environmental Lab(PMEL),
*  hereafter referred to as NOAA/PMEL/TMAP.
*
*  Access and use of this software shall impose the following
*  obligations and understandings on the user. The user is granted the
*  right, without any fee or cost, to use, copy, modify, alter, enhance
*  and distribute this software, and any derivative works thereof, and
*  its supporting documentation for any purpose whatsoever, provided
*  that this entire notice appears in all copies of the software,
*  derivative works and supporting documentation.  Further, the user
*  agrees to credit NOAA/PMEL/TMAP in any publications that result from
*  the use of this software or in any product that includes this
*  software. The names TMAP, NOAA and/or PMEL, however, may not be used
*  in any advertising or publicity to endorse or promote any products
*  or commercial entity unless specific written permission is obtained
*  from NOAA/PMEL/TMAP. The user also understands that NOAA/PMEL/TMAP
*  is not obligated to provide the user with any support, consulting,
*  training or assistance of any kind with regard to the use, operation
*  and performance of this software nor to provide the user with any
*  updates, revisions, new versions or "bug fixes".
*
*  THIS SOFTWARE IS PROVIDED BY NOAA/PMEL/TMAP "AS IS" AND ANY EXPRESS
*  OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
*  ARE DISCLAIMED. IN NO EVENT SHALL NOAA/PMEL/TMAP BE LIABLE FOR ANY
*  SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
*  RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
*  CONTRACT, NEGLIGENCE OR OTHER TORTUOUS ACTION, ARISING OUT OF OR IN
*  CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
*
*  Define an average lag axis using the provided calendar F and T axes.
*
*  The lag times in T for a time-of-forecast F0 of a forecast (dataset)
*  is defined as T - F0 where T are the forecasted times made in that
*  forecast.  The average lag T axis is then defined as an average of
*  these lag times in T for each time of forecast.
*
*  The lag times in F for a forecasted time T0 is defined as T0 - F
*  where F are the times-of-forecast for forecasts (datasets) that
*  forecasted that T0 time.  The average lag F axis is then defined
*  as an average of these lag times in F for each forecasted time.
*
*  Note that only one of these lag axes is created by this method;
*  which one is determined by the presences of '@LAG' in fmatch or
*  tmatch.  The return axes axlagf and axlagt are initialized to
*  unspecified_int4.  On successful return, one of those will be
*  a valid axis number.
*
*
* calling arguments
*   nsets          - number of datasets given to be aggregated (input)
*   member_sets    - array of open dataset numbers (input)
*   memb_filenames - array of dataset filenames (input)
*   fmatch         - method of matching times to fline (input, can be empty);
*                    eg, '@ASN', '@BIN', '@XACT', '@XACT:0.001', '@LAG', '@ASN@LAG', ...
*   findices       - index of each dataset on the time-of-forecast (F) axis (input)
*   tmatch         - method of matching times to tline (input, can be empty);
*                    eg, '@ASN', '@BIN', '@XACT', '@XACT:0.001', '@LAG', '@ASN@LAG', ...
*   t1indices      - index of each first-forecasted-time on the forecasted-time (T) axis (input)
*   tdims          - number of forecasted times for each forecast dataset (input)
*   axcalf         - time-of-forecast (F) axis used (input)
*   axcalt         - forecasted-times (T) axis used (input)
*   axlagf         - time-of-forecast lag (F) axis (output)
*   axlagt         - forecasted-times lag (T) axis (output)
*   status         - return status (output)

        IMPLICIT NONE

#include "tmap_dset.parm"
        include 'tmap_dims.parm'
        include 'ferret.parm'
        include 'errmsg.parm'
        include 'xdset_info.cmn_text'
        external xdset_info_data
        include 'xdset_info.cd_equiv'
        include 'xtm_grid.cmn_text'
        include 'xunits.cmn_text'
        include 'xvariables.cmn'
        include 'netcdf.inc'

* Calling argument declarations:
        INTEGER nsets, member_sets(nsets), findices(nsets),
     .          t1indices(nsets), tdims(nsets),
     .          axcalf, axcalt, axlagf, axlagt, status
        CHARACTER*512 memb_filenames(nsets)
        CHARACTER*(*) fmatch, tmatch

* Local variables:
        REAL*8  GET_LINE_COORD, TM_WORLD_AX, SECS_FROM_BC,
     .          tcoord, fcoord, lag, avlag, lastavlag, delta, first,
     .          lobnd, hibnd, udst, depsilon, epsil, bnd_epsil
        INTEGER*8 wksiz
        INTEGER TM_FIND_LIKE_DYN_LINE,
     .          fsiz, fsizp, tsiz, tsizp, wkspc, fidx, tidx, idx,
     .          sizi, sizip, sizj, units, lagax, idxi, idxj, nlags,
     .          pre_exist
        LOGICAL make_lagf, make_lagt
        CHARACTER nametmplt*8, linedir*2, calname*32, buff*256

        PARAMETER( depsilon = 2.0**(-42) )

        axlagf = unspecified_int4
        axlagt = unspecified_int4

        make_lagf = ( INDEX(fmatch, '@LAG') .GT. 0 )
        make_lagt = ( INDEX(fmatch, '@LAG') .GT. 0 )
*       The follow error should already have been caught
        IF ( make_lagf .AND. make_lagt )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .           'lag axes requested for both the F and T axes', *5000)

        wkspc = unspecified_int4
        fsiz  = line_dim(axcalf)
        fsizp = fsiz + 1
        tsiz  = line_dim(axcalt)
        tsizp = tsiz + 1
*       Two blocks needed; fsiz * tsiz for coordinates, (fsiz+1)*(tsiz+1) for bounds (actually excess)
        wksiz = fsiz * tsiz + fsizp * tsizp
        CALL CREATE_DYN_WRKSPC(wksiz, wkspc, status)
        IF ( status .NE. ferr_ok )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .                 'out of memory for computing lag axis', *5000)

        IF ( make_lagt ) THEN
*          Compute T-lag coordinates and bounds in units of SECS_FROM_BC
*          Assign the T coordinates and bounds
           DO tidx = 1,tsiz
              tcoord = TM_WORLD_AX(tidx, axcalt, BOX_MIDDLE)
              DO fidx = 1,fsiz
                 idx = (fidx - 1) * tsiz + tidx
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, tcoord)
              ENDDO
              tcoord = TM_WORLD_AX(tidx, axcalt, BOX_LO_LIM)
              DO fidx = 1,fsiz
                 idx = fsiz * tsiz + (fidx - 1) * tsizp + tidx
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, tcoord)
              ENDDO
           ENDDO
*          Add the upper edge of the T axis
           tcoord = TM_WORLD_AX(tsiz, axcalt, BOX_HI_LIM)
           DO fidx = 1,fsiz
              idx = fsiz * tsiz + (fidx - 1) * tsizp + tsizp
              CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, tcoord)
           ENDDO
*          Subtract the F coordinates from the T coordinates and bounds
           DO fidx = 1,fsiz
              fcoord = TM_WORLD_AX(fidx, axcalf, BOX_MIDDLE)
              DO tidx = 1,tsiz
                 idx = (fidx - 1) * tsiz + tidx
                 tcoord = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
                 lag = tcoord - fcoord
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, lag)
              ENDDO
              DO tidx = 1,tsizp
                 idx = fsiz * tsiz + (fidx - 1) * tsizp + tidx
                 tcoord = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
                 lag = tcoord - fcoord
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, lag)
              ENDDO
           ENDDO
           sizi  = tsiz
           sizip = tsizp
           sizj  = fsiz
           nametmplt = 'TF_LAG_T'
           linedir = 'TT'
           calname = line_cal_name(axcalt)
           units = line_unit_code(axcalt)
        ELSE
*          Compute F-lag coordinates and bounds in units of SECS_FROM_BC
*          Note that each block of values are put into the workspace transpose to above
*          Assign the F coordinates and bounds
           DO fidx = 1,fsiz
              fcoord = TM_WORLD_AX(fidx, axcalf, BOX_MIDDLE)
              DO tidx = 1,tsiz
                 idx = (tidx - 1) * fsiz + fidx
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, fcoord)
              ENDDO
              fcoord = TM_WORLD_AX(fidx, axcalf, BOX_LO_LIM)
              DO tidx = 1,tsiz
                 idx = tsiz * fsiz + (tidx - 1) * fsizp + fidx
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, fcoord)
              ENDDO
           ENDDO
*          Add the upper edge of the F axis
           fcoord = TM_WORLD_AX(fsiz, axcalf, BOX_HI_LIM)
           DO tidx = 1,tsiz
              idx = tsiz * fsiz + (tidx - 1) * fsizp + fsizp
              CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, fcoord)
           ENDDO
*          Subtract the F coordinates and bounds from the T coordinates
           DO tidx = 1,tsiz
              tcoord = TM_WORLD_AX(tidx, axcalt, BOX_MIDDLE)
              DO fidx = 1,fsiz
                 idx = (tidx - 1) * fsiz + fidx
                 fcoord = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
                 lag = tcoord - fcoord
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, lag)
              ENDDO
              DO fidx = 1,fsizp
                 idx = tsiz * fsiz + (tidx - 1) * fsizp + fidx
                 fcoord = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
                 lag = tcoord - fcoord
                 CALL PUT_LINE_COORD(workmem(wkspc)%ptr, idx, lag)
              ENDDO
           ENDDO
           sizi  = fsiz
           sizip = fsizp
           sizj  = tsiz
           nametmplt = 'TF_LAG_F'
           linedir = 'FF'
           calname = line_cal_name(axcalf)
           units = line_unit_code(axcalf)
        ENDIF

*       Create the average lag axis; being time differences, an axis
*       of time units, but not actually a calendar axis - no T0
        CALL TM_ALLO_TMP_LINE(lagax, status)
        IF ( status .NE. ferr_ok )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .                 'unable to allocate a lag axis', *5300)
        CALL TM_NEW_LINE_NAME(nametmplt, buff)
        line_name(lagax)         = buff
        line_name_orig(lagax)    = buff
        line_dattype(lagax)      = nf_double
        line_modulo(lagax)       = .FALSE.
        line_modulo_len(lagax)   = 0.0
        line_shift_origin(lagax) = .FALSE.
        line_dim_only(lagax)     = .FALSE.
        line_cal_name(lagax)     = calname
        line_unit_code(lagax)    = units
        line_units(lagax)        = un_name(units)
        line_tunit(lagax)        = un_convert(units)
        line_direction(lagax)    = linedir
        line_t0(lagax)           = char_init20
*       start with irregular axis, then check later if regular
        line_regular(lagax)      = .FALSE.
        CALL GET_LINE_DYNMEM(nsets, lagax, status)
        IF ( status .NE. ferr_ok )
     .     CALL ERRMSG(ferr_aggregate_error, status,
     .                 'out of memory for average lag axis', *5300)

* TODO: the following is WRONG - the length of the lag axis is dictated by the longest 
*       T or F in each dataset; should be dictated by findices/nsets and t1indices/tdims
        line_dim(lagax)          = sizi

*       Compute the average lag axis.
*       Average the axis with coordinates I=firsti:lasti for each axis J in jindices
*       where idx = (J-1)*sizi + I is the index of the coordinate in the workspace.
*       Also average the edges I=firsti:lasti+1 for each axis J in jindices where
*       idx = sizi * sizj + (J-1)*sizip + I is the index of the edge in the workspace.

* TODO: the following blocks are WRONG - axes need to be offset and truncated to obtain 
*       matching lags times; should be dictated by findices/nsets and t1indices/tdims
*       Compute the lowest edge to initialize lastavlag
*       (the lowest lag edge could be negative)
        idxi = 1
        avlag = 0.0
        DO idxj = 1,sizj
           idx = sizi * sizj + (idxj - 1) * sizip + idxi
           lag = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
           avlag = avlag + lag
        ENDDO
        avlag = avlag / DBLE(sizj)
        CALL PUT_LINE_COORD(lineedg(lagax)%ptr, idxi, avlag)
        lastavlag = avlag

        DO WHILE ( idxi .LE. sizi )
*          Compute the coordinate
           avlag = 0.0
           DO idxj = 1,sizj
              idx = (idxj - 1) * sizi + idxi
              lag = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
              avlag = avlag + lag
           ENDDO
           avlag = avlag / DBLE(sizj)
*          Make sure lag coordinates and edges are monotonically increasing
           IF ( avlag .LE. lastavlag ) GOTO 5010
           CALL PUT_LINE_COORD(linemem(lagax)%ptr, idxi, avlag)
           lastavlag = avlag

           idxi = idxi + 1

*          Compute the upper edge
           avlag = 0.0
           DO idxj = 1,sizj
              idx = sizi * sizj + (idxj - 1) * sizip + idxi
              lag = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
              avlag = avlag + lag
           ENDDO
           avlag = avlag / DBLE(sizj)
*          Make sure lag edges/coordinates are monotonically increasing
           IF ( avlag .LE. lastavlag ) GOTO 5010
           CALL PUT_LINE_COORD(lineedg(lagax)%ptr, idxi, avlag)
           lastavlag = avlag
        ENDDO

* TODO: the following is WRONG - axes need to be offset and truncated as above
*       Verify the coordinates of each individual lag axis fall
*       in the appropriate average lag axis cell
        lobnd = GET_LINE_COORD(lineedg(lagax)%ptr, 1)
        DO idxi = 1,sizi
           hibnd = GET_LINE_COORD(lineedg(lagax)%ptr, idxi+1)
           DO idxj = 1,sizj
              idx = (idxj - 1) * sizi + idxi
              lag = GET_LINE_COORD(workmem(wkspc)%ptr, idx)
              IF ( (lag .LT. lobnd) .OR. (lag .GT. hibnd) ) GOTO 5010
           ENDDO
           lobnd = hibnd
        ENDDO


* The rest of this should be correct
        nlags = line_dim(lagax)

*       Convert into units of this lag axis
        udst = un_convert(units)
        DO idx = 1,nlags
           avlag = GET_LINE_COORD(linemem(lagax)%ptr, idx)
           avlag = avlag / udst
           CALL PUT_LINE_COORD(linemem(lagax)%ptr, idx, avlag)
        ENDDO
        DO idx = 1,nlags+1
           avlag = GET_LINE_COORD(lineedg(lagax)%ptr, idx)
           avlag = avlag / udst
           CALL PUT_LINE_COORD(lineedg(lagax)%ptr, idx, avlag)
        ENDDO

*       Check if this is a regular axis.
*       Since in units of the axis, use standard small error value
        epsil = depsilon
*       A little more leniency for edges
        bnd_epsil = 16.0 * depsilon
        first = GET_LINE_COORD(linemem(lagax)%ptr, 1)
        IF ( nlags .GT. 1 ) THEN
           lag = GET_LINE_COORD(linemem(lagax)%ptr, nlags)
           delta = (lag - first) / (nlags - 1.0)
           IF ( nlags .GT. 2 ) THEN
              DO idx = 2,nlags-1
                 lag = GET_LINE_COORD(linemem(lagax)%ptr, idx)
                 IF ( ABS(first + (idx - 1.0) * delta - lag) 
     .                    .GT. epsil ) GOTO 100
              ENDDO
           ENDIF
*          Check if the provided cell bounds are at the "midpoints".
           DO idx = 1,nlags+1
              lag = GET_LINE_COORD(lineedg(lagax)%ptr, idx)
              IF ( ABS(first + (idx - 1.5) * delta - lag) 
     .                 .GT. bnd_epsil ) GOTO 100
           ENDDO
        ELSE
           lobnd = GET_LINE_COORD(lineedg(lagax)%ptr, 1)
           hibnd = GET_LINE_COORD(lineedg(lagax)%ptr, 2)
           IF ( ABS(first - 0.5 * (lobnd + hibnd)) 
     .              .GT. bnd_epsil ) GOTO 100
           delta = hibnd - lobnd
        ENDIF

*       If we got here, the axis is regular (including edges)
        CALL FREE_LINE_DYNMEM(lagax)
        line_regular(lagax) = .TRUE.
        line_start(lagax) = first
        line_delta(lagax) = delta

  100   pre_exist = TM_FIND_LIKE_DYN_LINE(lagax)
        IF ( pre_exist .NE. unspecified_int4 ) THEN
*          Axis is exact duplicate of existing axis; use that axis
           IF ( .NOT. line_regular(lagax) )
     .           CALL FREE_LINE_DYNMEM(lagax)
           lagax = pre_exist
        ELSE
*          Move the tmp_line to a dyn_line
           CALL TM_RE_ALLO_TMP_LINE(lagax)
        ENDIF

*       Success
        IF ( make_lagt ) THEN
           axlagt = lagax
        ELSE
           axlagf = lagax
        ENDIF
 1000   status = ferr_ok
 5000   RETURN

 5010   CALL ERRMSG(ferr_aggregate_error, status,
     .              'date axis is too irregular to generate '//
     .              'a reasonable average lag axis', *5200)

 5200   CALL FREE_LINE_DYNMEM(lagax)
 5300   IF ( wkspc .NE. unspecified_int4 )
     .     CALL FREE_WS_DYNMEM(wkspc)
        RETURN

        END

