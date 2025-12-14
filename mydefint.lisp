(defun classify-imag-part (z)
 "Classifies the expression `z` based on the sign of its imaginary part.

Returns:
  'lower     — if the imaginary part is negative
  'realaxis  — if the imaginary part is zero
  'upper     — if the imaginary part is positive

Uses `asksign` to determine the sign."
  ;; It's unlikely that `asksign will throw an error, but it's possible, for example `asksign(und)`. So
  ;; use errcatch.
  (case (car (errcatch ($asksign ($imagpart z))))
    ($neg   'lower)
    ($zero  'realaxis)
    ($pos   'upper)
    (otherwise (merror (intl:gettext "classify-imag-part: Unexpected result from asksign")))))

(defun defint-by-symmetric-interval (e x lo hi)
"When the interval is symmetric and the integrand odd, return 0; otherwise, return nil."
    (if (and (or (and (eq lo '$minf) (eq hi '$inf))
	               (alike1 hi (mul -1 lo)))
		         (oddfn e x)) 
		0
		nil))

(defun defint-by-mbag (e x lo hi)
"If e is an mbag, map q -> integrate(q,x,lo,hi) over e; if not, return nil."
	(and (mbagp e)
	     (fapply (caar e) (mapcar #'(lambda (q) ($mydefint q x lo hi)) (cdr e)))))

(defun defint-rational-minf-to-inf (w x lo hi)
"Integrates the rational function x -> w over the real line from minf to inf using Cauchy's integral
 formula. If a pole lies on the real axis, print a messages that the result is a principal value 
 integral. When x -> w isn't a rational function, return nil."
  (let ((pv nil) ($algebraic t) (cw 0) (sum nil) (p ($num w)) (q ($denom w)))
    (cond
      ((not (ratp w x)) nil)
      ((or (not (eq lo '$minf)) (not (eq hi '$inf))) nil)
      ((eq t (mgrp (add 1 ($hipow ($expand p) x)) ($hipow ($expand q) x)))
       (diverg))  ; Degree condition: integral diverges
      (t
       (multiple-value-bind (poles mp) (solve-with-multiplicities q x)
         (flet ((residue-weight (p)
                 (let ((plane (classify-imag-part p)))
                   (cond ((eq plane 'lower) 0)
                         ((eq plane 'realaxis)
                          (setq pv t)
                          (div 1 2))
                         (t 1)))))
		   ;; weighted sum of residues
       (dolist (pk poles)
          (setq cw (residue-weight pk))
          (if (not (eql 0 cw))
            (push (mul cw (residue-by-methods w x pk)) sum)))
		  ;; When a pole is on the x-axis, inform the alert reader that this is a principal value integral
      (when pv
          (mtell (intl:gettext "Principal value ~%")))
      ($expand (mul 2 '$%i '$%pi (fapply 'mplus sum)) 0 0)))))))

(defun defint-rational-zero-to-inf (w x lo hi)
"Integrates the rational function x -> w over the real line from minf to inf using Cauchy's integral
 formula. If a pole lies on the real axis, print a messages that the result is a principal value 
 integral. When x -> w isn't a rational function, return nil."
  (let ((pv nil) ($algebraic t) (cw 0) (sum nil) (p ($num w)) (q ($denom w)))
    (cond
      ((not (ratp w x)) nil)
      ((or (not (eql lo 0)) (not (eq hi '$inf))) nil)
      ((eq t (mgrp (add 1 ($hipow ($expand p) x)) ($hipow ($expand q) x)))
       (diverg))  ; Degree condition: integral diverges
      (t
       (setq w (mul (ftake '%plog (neg x)) w))
       (multiple-value-bind (poles mp) (solve-with-multiplicities q x)
         (mtell "poles = ~M ~%" (fapply'mlist poles))
         (flet ((residue-weight (p)
                 (let ((plane (classify-imag-part p)))
                   (cond ((eq plane 'lower) 1)
                         ((eq plane 'realaxis)
                          (setq pv t)
                          (div 1 2))
                         (t 1)))))
		   ;; weighted sum of residues
       (dolist (pk poles)
          (setq cw (residue-weight pk))
          (mtell "cw = ~M ~%" cw)
          (if (not (eql 0 cw))
            (push (mul cw (residue-by-methods w x pk)) sum)))
		  ;; When a pole is on the x-axis, inform the alert reader that this is a principal value integral
      (when pv
          (mtell (intl:gettext "Principal value ~%")))

      ($expand (mul -1 (fapply 'mplus sum)) 1 0)))))))

(defun defint-by-nounform (e x lo hi) 
   (setq ans (ftake '%integrate e x lo hi)))

(defvar *defint-method-info* t)
(defvar *defint-methods* 
      (list 'defint-by-mbag
            'defint-by-symmetric-interval 
            'defint-rational-minf-to-inf
            'defint-rational-zero-to-inf
            'defint-by-nounform))

(defun $mydefint (e x lo hi)
 (catch 'finished
    (dolist (fn *defint-methods*)
      (let ((ans (funcall fn e x lo hi)))
        (when ans
          (when *defint-method-info*
            (mtell (intl:gettext "Definite integration method: ~A succeeded.~%") fn))
          (throw 'finished ans))))))

(defun defint-xxx (e x lo hi)
  (let (($logabs nil)
        ($exptsubst t) ;not sure why?
        (limitp t)
        ($noprincipal nil)
        ($expintrep '$gamma_incomplete) ;not sure why?
        (xx ($gensym "z"))
		    (ans)
		    ($errormsg nil)
        (cntx ($supcontext)))

    ;; We don't want asksign questions about a gensym, so declare xx to be internal
    (putprop xx t 'internal)

    (when (eq t (mgrp lo hi))
      (setq e ($limit (neg e)))
      (psetf lo hi hi lo))

    (unwind-protect
        (progn
          (setq e (maxima-substitute xx x e))
          (assume (ftake 'mgreaterp xx lo))
          (assume (ftake 'mgreaterp hi xx))
          (assume (ftake 'mgreaterp hi lo))
          ;; Vital resimplification, e.g. for integrate(x*signum(x^3 - 8), x, 2, 5)
          (setq e (resimplify e))
          (setq e (%i-out-of-denom e)) ;not sure about this
          ;(setq exp (tansc-var exp xx)) ;not sure about this either
		      (setq e ($gfactor e))
          (catch 'finished ans
		          (dolist (fn *defint-methods*)
		    	       (setq ans (funcall fn e xx lo hi))
			           (when ans 
			              (throw 'finished ans))))
          (when (not (freeof '$ind ans))
		   	      (setq ans (ftake '%integrate e x lo hi)))

		      (when (not (free ans xx))
		  	     (setq ans (maxima-substitute x xx ans)))
 
          ans)
      ($killcontext cntx))))