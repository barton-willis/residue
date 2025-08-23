;;; -*-  Mode: Lisp; Package: Maxima; Syntax: Common-Lisp; Base: 10 -*- ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;     The data in this file contains enhancements.                   ;;;;;
;;;                                                                    ;;;;;
;;;  Copyright (c) 1984,1987 by William Schelter,University of Texas   ;;;;;
;;;     All rights reserved                                            ;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;     (c) Copyright 1982 Massachusetts Institute of Technology         ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :maxima)

(macsyma-module residu)

(load-macsyma-macros rzmac)

(declare-top (special $noprincipal 
          $algexact $solveexplicit $solvedecomposes
		      var *roots *failures
		      sn* sd* genvar zn))


;; Compute the poles (roots) of the polynomial D and return them.
;; Divide these roots into several parts: Those in REGION, REGION1,
;; and everywhere else.  These are returned in a list.  (In a more
;; modern style, we'd probably return them in 4 different values.)
;;
;; The regions are determined by functions REGION and REGION1, which
;; should return non-NIL if the root is in the given region.
;;
;; The description below applies if *SEMIRAT* is NIL.  If *SEMIRAT* is
;; non-NIL, somewhat different results are returned.  I (rtoy) am not
;; exactly sure what *SEMIRAT* is intended to mean.
;;
;; The first part of the list of the form ((r1 (x - r1)^d1) (r2 (x -
;; r2)^d2) ...) where r1, r2 are the roots, d1, d2 are the
;; multiplicity of each root, and x is the variable.
;;
;; The second part is a list of the repeated roots in REGION.  Each
;; element of the list is of the form (r d) where r is the root and d
;; is the multiplicity.
;;
;; The third part is a list of the simple roots in REGION.
;;
;; Finally, the fourth part is NIL, unless *semirat* is T.
(defun polelist (d region region1)
  (prog (roots $breakup r rr ss r1 s pole wflag cf)
     (setq wflag t)
     (setq *leadcoef* (polyinx d var 'leadcoef))
     (setq roots (solvecase d))
     (if (eq roots 'failure) (return ()))
     ;; Loop over all the roots.  SOLVECASE returns the roots in the form
     ;; ((x = r1) mult1
     ;;  (x = r2) mult2
     ;;  ...)

   loop1
     (cond ((null roots)
	    (cond ((and *semirat* nil
			(> (+ (length s) (length r))
			   (+ (length ss) (length rr))))
		   ;; Return CF, repeated roots (*semirat*), simple
		   ;; roots (*semirat*), roots in region 1.
		   (return (list cf rr ss r1)))
		  (t
		   ;; Return CF, repeated roots, simple roots, roots in region 1.
		   (return (list cf r s r1)))))
	   (t
	    ;; Set POLE to the actual root and set D to the
	    ;; multiplicity of the root.
	    (setq pole (caddar roots))
	    (setq d (cadr roots))
	    (cond (*leadcoef*
		   ;; Is it possible for LEADCOEF to be NIL ever?
		   ;;
		   ;; Push (pole (x - pole)^d) onto the list CF.
		   (setq cf (cons pole
				  (cons
				   (m^ (m+ var (m* -1 pole))
				       d)
				   cf)))))))
     ;; Don't change the order of clauses here.  We want to call REGION and then REGION1.
     (cond ((funcall region pole)
	    ;; The pole is in REGION
	    (cond ((equal d 1)
		   ;; A simple pole, so just push the pole onto the list S.
		   (push pole s))
		  (t
		   ;; A multiple pole, so push (pole d) onto the list R.
		   (push (list pole d) r))))
	   ((funcall region1 pole)
	    ;; The pole is in REGION1
	    (cond ((not $noprincipal)
		   ;; Put the pole onto the R1 list.  (Don't know what
		   ;; $NOPRINCIPAL is.)
		   (push pole r1))
		  (t
		   ;; Return NIL if we get here.
		   (return nil))))
	   (*semirat*
	    ;; (What does *SEMIRAT* mean?)  Anyway if we're here, the
	    ;; pole is not in REGION or REGION1, so push the pole onto
	    ;; SS or RR depending if the pole is repeated or not.
	    (cond ((equal d 1)
		   (push pole ss))
		  (t (push (list pole d) rr)))))
     ;; Pop this root and multiplicity and move on.
     (setq roots (cddr roots))
     (go loop1)))

;; Like POLELIST, but dependency on VAR is explicit.  Use this instead
;; when possible.
(defun polelist-var (var1 d region region1)
  (prog (roots $breakup r rr ss r1 s pole wflag cf)
     (setq wflag t)
     (setq *leadcoef* (polyinx d var1 'leadcoef))
     (setq roots (solvecase-var d var1))
     (if (eq roots 'failure) (return ()))
     ;; Loop over all the roots.  SOLVECASE returns the roots in the form
     ;; ((x = r1) mult1
     ;;  (x = r2) mult2
     ;;  ...)

   loop1
     (cond ((null roots)
	    (cond ((and *semirat*
			(> (+ (length s) (length r))
			   (+ (length ss) (length rr))))
		   ;; Return CF, repeated roots (*semirat*), simple
		   ;; roots (*semirat*), roots in region 1.
		   (return (list cf rr ss r1)))
		  (t
		   ;; Return CF, repeated roots, simple roots, roots in region 1.
		   (return (list cf r s r1)))))
	   (t
	    ;; Set POLE to the actual root and set D to the
	    ;; multiplicity of the root.
	    (setq pole (caddar roots))
	    (setq d (cadr roots))
	    (cond (*leadcoef*
		   ;; Is it possible for LEADCOEF to be NIL ever?
		   ;;
		   ;; Push (pole (x - pole)^d) onto the list CF.
		   (setq cf (cons pole
				  (cons
				   (m^ (m+ var1 (m* -1 pole))
				       d)
				   cf)))))))
     ;; Don't change the order of clauses here.  We want to call REGION and then REGION1.
     (cond ((funcall region pole)
	    ;; The pole is in REGION
	    (cond ((equal d 1)
		   ;; A simple pole, so just push the pole onto the list S.
		   (push pole s))
		  (t
		   ;; A multiple pole, so push (pole d) onto the list R.
		   (push (list pole d) r))))
	   ((funcall region1 pole)
	    ;; The pole is in REGION1
	    (cond ((not $noprincipal)
		   ;; Put the pole onto the R1 list.  (Don't know what
		   ;; $NOPRINCIPAL is.)
		   (push pole r1))
		  (t
		   ;; Return NIL if we get here.
		   (return nil))))
	   (*semirat*
	    ;; (What does *SEMIRAT* mean?)  Anyway if we're here, the
	    ;; pole is not in REGION or REGION1, so push the pole onto
	    ;; SS or RR depending if the pole is repeated or not.
	    (cond ((equal d 1)
		   (push pole ss))
		  (t (push (list pole d) rr)))))
     ;; Pop this root and multiplicity and move on.
     (setq roots (cddr roots))
     (go loop1)))

(defun solvecase (e)
  (cond ((not (among var e)) nil)
	(t (let (*failures *roots)
	     (solve e var 1)
	     (cond (*failures 'failure)
		   ((null *roots) ())
		   (t *roots))))))

(defun solvecase-var (e var1)
  (cond ((not (among var1 e)) nil)
	(t (let (*failures *roots)
	     (solve e var1 1)
	     (cond (*failures 'failure)
		   ((null *roots) ())
		   (t *roots))))))

;; Compute the sum of the residues of n/d.
(defun res (n d region region1)
  (let ((pl (polelist d region region1))
	dp a b c factors *leadcoef*)
    (cond
      ((null pl) nil)
      (t
       (setq factors (car pl))
       (setq pl (cdr pl))
       ;; PL now contains the list of the roots in region, roots in
       ;; region1, and everything else.
       (cond ((or (cadr pl)
		  (caddr pl))
	      (setq dp (sdiff d var))))
       (cond ((car pl)
	      ;; Compute the sum of the residues of n/d for the
	      ;; multiple roots in REGION.
	      (setq a (m+l (residue n (cond (*leadcoef* factors)
					    (t d))
				    (car pl)))))
	     (t (setq a 0)))
       (cond ((cadr pl)
	      ;; Compute the sum of the residues of n/d for the simple
	      ;; roots in REGION1.  Since the roots are simple, we can
	      ;; use RES1 to compute the residues instead of the more
	      ;; complicated $RESIDUE.  (This works around a bug
	      ;; described in bug 1073338.)
	      #+nil
	      (setq b (m+l (mapcar #'(lambda (pole)
				       ($residue (m// n d) var pole))
				   (cadr pl))))
	      (setq b (m+l (res1 n dp (cadr pl)))))
	     (t (setq b 0)))
       (cond ((caddr pl)
	      ;; Compute the sum of the residues of n/d for the roots
	      ;; not in REGION nor REGION1.
	      (setq c (m+l (mapcar #'(lambda (pole)
				       ($residue (m// n d) var pole))
				   (caddr pl)))))
	     (t (setq c ())))
       ;; Return the sum of the residues in the two regions and the
       ;; sum of the residues outside the two regions.
       (list (m+ a b) c)))))

(defun res-var (var1 n d region region1)
  (let ((pl (polelist-var var1 d region region1))
	dp a b c factors *leadcoef*)
    (cond
      ((null pl) nil)
      (t
       (setq factors (car pl))
       (setq pl (cdr pl))
       ;; PL now contains the list of the roots in region, roots in
       ;; region1, and everything else.
       (cond ((or (cadr pl)
		  (caddr pl))
	      (setq dp (sdiff d var1))))
       (cond ((car pl)
	      ;; Compute the sum of the residues of n/d for the
	      ;; multiple roots in REGION.
	      (setq a (m+l (let ((var var1))
                             ;; This binding is very important for
                             ;; defint.  It prevents the plog
                             ;; simplifier from simplifying plog until
                             ;; the pole has been substituted in when
                             ;; computing the residue.
                             (declare (special var))
                             (residue-var var1 n (cond (*leadcoef* factors)
					               (t d))
				          (car pl))))))
	     (t (setq a 0)))
       (cond ((cadr pl)
	      ;; Compute the sum of the residues of n/d for the simple
	      ;; roots in REGION1.  Since the roots are simple, we can
	      ;; use RES1 to compute the residues instead of the more
	      ;; complicated $RESIDUE.  (This works around a bug
	      ;; described in bug 1073338.)
	      #+nil
	      (setq b (m+l (mapcar #'(lambda (pole)
				       ($residue (m// n d) var1 pole))
				   (cadr pl))))
	      (setq b (m+l (res1-var var1 n dp (cadr pl)))))
	     (t (setq b 0)))
       (cond ((caddr pl)
	      ;; Compute the sum of the residues of n/d for the roots
	      ;; not in REGION nor REGION1.
	      (setq c (m+l (mapcar #'(lambda (pole)
				       ($residue (m// n d) var1 pole))
				   (caddr pl)))))
	     (t (setq c ())))
       ;; Return the sum of the residues in the two regions and the
       ;; sum of the residues outside the two regions.
       (list (m+ a b) c)))))

(defun residue (zn factors pl)
  (residue-var var zn factors pl))

(defun residue-var (var1 zn factors pl)
  (cond (*leadcoef*
	 (mapcar #'(lambda (j)
		     (destructuring-let (((factor1 factor2)
                                          (remfactor factors (car j) zn)))
		       (resm0-var var1 factor1 factor2 (car j) (cadr j))))
		 pl))
	(t (mapcar #'(lambda (j)
		       (resm1-var var1 (div zn factors) (car j)))
		   pl))))

;; Compute the residues of zn/d for the simple poles in the list PL1.
;; The poles must be simple because ZD must be the derivative of
;; denominator.  For simple poles the residue can be computed as
;; limit(n(z)/d(z)*(z-a),z,a).  Since the pole is simple we have the
;; indeterminate form (z-a)/d(z).  Use L'Hospital's rule to get
;; 1/d'(z).  Hence, the residue is easily computed as n(a)/d'(a).
(defun res1 (zn zd pl1)
  (setq zd (div* zn zd))
  (mapcar #'(lambda (j)
	      ;; In case the pole is messy, call $RECTFORM.  This
	      ;; works around some issues with gcd bugs in certain
	      ;; cases.  (See bug 1073338.)
        ;(subin j zd))
        (maxima-substitute j var zd))
	      ;($rectform ($expand (subin ($rectform j) zd))))
	  pl1))

(defun res1-var (var1 zn zd pl1)
  (let ((w (div zn zd)))
     (mapcar #'(lambda (p) (maxima-substitute p var1 w)) pl1)))

(defun resprog0 (f g n n2)
  (prog (a b c r)
     (setq a (resprog f g))
     (setq b (cadr a) c (ptimes (cddr a) n2) a (caar a))
     (setq a (ptimes n a) b (ptimes n b))
     (setq r (pdivide a g))
     (setq a (cadr r) r (car r))
     (setq b (cons (pplus (ptimes (car r) f) (ptimes (cdr r) b))
		   (cdr r)))
     (return (cons (cons (car a) (ptimes (cdr a) c))
		   (cons (car b) (ptimes (cdr b) c))))))

(defun resprog0-var (var1 f g n n2)
  (prog (a b c r)
     (setq a (resprog-var var1 f g))
     (setq b (cadr a) c (ptimes (cddr a) n2) a (caar a))
     (setq a (ptimes n a) b (ptimes n b))
     (setq r (pdivide a g))
     (setq a (cadr r) r (car r))
     (setq b (cons (pplus (ptimes (car r) f) (ptimes (cdr r) b))
		   (cdr r)))
     (return (cons (cons (car a) (ptimes (cdr a) c))
		   (cons (car b) (ptimes (cdr b) c))))))


(defun resm0 (e n pole m)
  (setq e (div n e))
  (setq e ($diff e var (1- m)))
  (setq e ($rectform ($expand (subin pole e))))
  (div e (simplify `((mfactorial) ,(1- m)))))

(defun resm0-var (var1 e n pole m)
  (setq e (div n e))
  (setq e ($diff e var1 (1- m)))
  (setq e ($rectform ($expand (subin-var pole e var1))))
  (div e (simplify `((mfactorial) ,(1- m)))))

(defun remfactor (l p n)
  (prog (g)
   loop (cond ((null l)
	       (return (list (m*l (cons *leadcoef* g)) n)))
	      ((equal p (car l)))
	      (t (setq g (cons (cadr l) g))))
   (setq l (cddr l))
   (go loop)))

(defun resprog (p1b p2b)
  (prog (temp coef1r coef2r fac coef1s coef2s zeropolb f1 f2)
     (setq coef2r (setq coef1s 0))
     (setq coef2s (setq coef1r 1))
     b1   (cond ((not (< (pdegree p1b var) (pdegree p2b var))) (go b2)))
     (setq temp p2b)
     (setq p2b p1b)
     (setq p1b temp)
     (setq temp coef2r)
     (setq coef2r coef1r)
     (setq coef1r temp)
     (setq temp coef2s)
     (setq coef2s coef1s)
     (setq coef1s temp)
     b2   (cond ((zerop (pdegree p2b var))
		 (return (cons (cons coef2r p2b) (cons coef2s p2b)))))
     (setq zeropolb (psimp var
			   (list (- (pdegree p1b var) (pdegree p2b var))
				 1)))
     (setq fac (pgcd (caddr p1b) (caddr p2b)))
     (setq f1 (pquotient (caddr p1b) fac))
     (setq f2 (pquotient (caddr p2b) fac))
     (setq p1b (pdifference (ptimes f2 (psimp (car p1b) (cdddr p1b)))
			    (ptimes f1
				    (ptimes zeropolb
					    (psimp (car p2b)
						   (cdddr p2b))))))
     (setq coef1r (pdifference (ptimes f2 coef1r)
			       (ptimes (ptimes f1 coef2r) zeropolb)))
     (setq coef1s (pdifference (ptimes f2 coef1s)
			       (ptimes (ptimes f1 coef2s) zeropolb)))
     (go b1)))

(defun resprog-var (var1 p1b p2b)
  (prog (temp coef1r coef2r fac coef1s coef2s zeropolb f1 f2)
     (setq coef2r (setq coef1s 0))
     (setq coef2s (setq coef1r 1))
     b1   (cond ((not (< (pdegree p1b var1) (pdegree p2b var1))) (go b2)))
     (setq temp p2b)
     (setq p2b p1b)
     (setq p1b temp)
     (setq temp coef2r)
     (setq coef2r coef1r)
     (setq coef1r temp)
     (setq temp coef2s)
     (setq coef2s coef1s)
     (setq coef1s temp)
     b2   (cond ((zerop (pdegree p2b var1))
		 (return (cons (cons coef2r p2b) (cons coef2s p2b)))))
     (setq zeropolb (psimp var1
			   (list (- (pdegree p1b var1) (pdegree p2b var1))
				 1)))
     (setq fac (pgcd (caddr p1b) (caddr p2b)))
     (setq f1 (pquotient (caddr p1b) fac))
     (setq f2 (pquotient (caddr p2b) fac))
     (setq p1b (pdifference (ptimes f2 (psimp (car p1b) (cdddr p1b)))
			    (ptimes f1
				    (ptimes zeropolb
					    (psimp (car p2b)
						   (cdddr p2b))))))
     (setq coef1r (pdifference (ptimes f2 coef1r)
			       (ptimes (ptimes f1 coef2r) zeropolb)))
     (setq coef1s (pdifference (ptimes f2 coef1s)
			       (ptimes (ptimes f1 coef2s) zeropolb)))
     (go b1)))

;;;Looks for polynomials. puts polys^(pos-num) in sn* polys^(neg-num) in sd*.
(defun snumden (e)
  (cond ((or (atom e)
	     (mnump e))
	 (setq sn* (cons e sn*)))
	((and (mexptp e)
	      (integerp (caddr e)))
	 (cond ((polyinx (cadr e) var nil)
		(cond ((minusp (caddr e))
		       (setq sd* (cons (cond ((equal (caddr e) -1) (cadr e))
					     (t (m^ (cadr e)
						    (- (caddr e)))))
				       sd*)))
		      (t (setq sn* (cons e sn*)))))))
	((polyinx e var nil)
	 (setq sn* (cons e sn*)))))

;; Like SNUMDEN, but dependency on VAR is explicit.  Use this instead
;; when possible.
(defun snumden-var (e var1 sn sd)
  (cond ((or (atom e)
	     (mnump e))
	 (setq sn (cons e sn))
         (values t sn sd))
	((and (mexptp e)
	      (integerp (caddr e)))
	 (cond ((polyinx (cadr e) var1 nil)
		(cond ((minusp (caddr e))
		       (setq sd (cons (cond ((equal (caddr e) -1) (cadr e))
					     (t (m^ (cadr e)
						    (- (caddr e)))))
				      sd))
                       (values t sn sd))
		      (t
                       (setq sn (cons e sn))
                       (values t sn sd))))))
	((polyinx e var1 nil)
	 (setq sn (cons e sn))
         (values t sn sd))))

(setq sn* nil sd* nil)

(defvar *residue-method-info* t
  "If non-nil, residue methods may print informational messages when they succeed.")

(defvar *residue-methods* nil)

;; The methods residue-nounform and residue-by-simp are required, but the other methods are 
;; optional.
(eval-when (:load-toplevel :execute)
  (setq *residue-methods*
        (list 'residue-by-misc-undefined 
              'residue-by-freeof 
              'residue-by-infinity-transform
              'residue-rational
              'residue-by-taylor
              'residue-by-taylor-asym
              'residue-by-powerseries 
              'residue-by-simp
              'residue-nounform))
    "List of methods for computing residues. Each method should accept three arguments 
  `(e x pt)` and return either a valid residue or nil if it cannot compute the result.")            

(defun solve-with-multiplicities (q x)
  "Return a list of the solutions to the equation `q=0` for `x`. Return both a CL list
 of the solutions and a CL list of the coresponding multiplicities. When solve fails, 
 return nil. The solutions are not in the form of an equation. "
  (let* (($breakup nil)
         ($programmode t)
         ($solvedecomposes t)
         ($solveexplicit t)
         ($algexact t)
         ($solvefactors t)
         (ans (mapcar #'$rhs (cdr ($solve ($factor q) x)))))
    ;; When solutions are missing, return (values nil nil)
    (if (eql ($hipow ($expand q) x)
             (reduce #'+ (cdr $multiplicities)))
        (values ans (cdr $multiplicities))
        (values nil nil))))        

(defun get-taylor-coeff (e x pt n &optional (vanish nil))
  "Return the nth coefficient of the Taylor expansion of expression `e`
   in variable `x` about point `pt`. The optional argument `vanish` is a Maxima
   expression known to vanish. When `vanish` is non-nil, the Taylor coefficient
   simplifier attempts to replace `vanish` with zero. If `taylor` fails to return
   a sum of integer powers, return NIL."
  (let ((silent-taylor-flag t)
        ($taylordepth 8)         ; Better than the default
        ($radexpand nil)
        ($maxtayorder t)
        ($demoivre nil)
        ($taylor_simplifier
         (if vanish
             #'(lambda (s) ($factor ($ratsubst 0 vanish s)))
             #'(lambda (s) ($factor (resimplify s)))))
        ($taylor_logexpand nil))
    (let ((ee (catch 'taylor-catch ($taylor e x pt n))))
      ;; Check if `ee` is a legitimate Taylor polynomial.
          (if (and 
               (not (eql ee 0))
               ($polynomialp ee (ftake 'mlist x) #'(lambda (s) (freeof x s)) #'integerp)
               (alike1 ($taylorinfo ee) (ftake 'mlist (ftake 'mlist x pt n))))
             ($ratdisrep (pscoeff1 ee x n))
             nil))))

(defun all-residues-rational-fun (p q x)
  "Compute all the residues of the rational function p/q at the poles of q using Heaviside's formula.
Arguments:
  p — numerator of the rational function (Maxima expression)
  q — denominator of the rational function (Maxima expression)
  x — variable with respect to which residues are computed
Returns:
  Two values:
    - A list of residues at each pole of q
    - A list of the corresponding poles
Notes:
  This function requires that solve can find all the roots of `q`; if it is not able to find all
  the roots, it returns nil. This function does not find the full partial fraction decomposition, 
  it only finds the first order pole terms."

  ;; set `poles` to a list of the zeros of `q` and set `ord` to a list of their multiplicities
  (multiple-value-bind (poles ord) (solve-with-multiplicities q x)
    ;; set `ppoles` to a list of gensyms, one for each pole,
    ;; set `qq` to a list of the factors of q, with the corresponding gensym substituted for the pole,
    ;; eventually `res` will be a list of the residues.
	(cond ((null poles) nil)
	      (t
    (let* ((res nil)
		   (qe ($expand q))
		   (lc (coeff qe x ($hipow qe x))) ;leading coefficient of q
		   (qq (mapcar #'(lambda (a b) (ftake 'mexpt (sub x a) b)) poles ord))
		   (qqq (mul lc (fapply 'mtimes qq))))
          (cond ((null poles) nil) ; return nil when unable to solve
              (t
               ;; For each factor qk of q, say qk = (x - pk)^nk, remove qk from q, call this qs
               ;; and evaluate diff(p/qs, nk-1) / (nk-1)! at x = pk. This is Heaviside’s formula.
			   (mapcar #'(lambda (qk pk)
			      (push (get-taylor-coeff
				          (div p (maxima-substitute 1 qk qqq)) 
						  x 
						  pk 
				          (if (mexptp qk) (sub (third qk) 1) 0)) res)) qq poles)
               (values (reverse res) poles))))))))

;; See https://math.stackexchange.com/questions/3480929/theorem-on-residue-of-composite-function
(defvar *residue* nil)
(defmfun $residue (e var p)
  (push (ftake 'mlist e var p) *residue*)
  ;; Error when var isn't a mapatom
  (when (not ($mapatom var))
    (merror (intl:gettext "residue: second argument must be a mapatom; found ~M") var))
  ;; Error when `p` depends on `var`
  (when (not (freeof var p))
    (merror (intl:gettext "residue: third argument must not depend on second argument; found ~M") p))
  ;; Error when `e` is an mbag
  (when (mbagp e)
	(merror (intl:gettext "residue: first argument must not be an mbag; found ~M") e))
   ;; Error when `p` is an mbag
  (when (mbagp p)
	(merror (intl:gettext "residue: third argument must not be an mbag; found ~M") p))

  ;; Unfortunately, when running the testsuite, sometimes `e` is not simplified. So, let us
  ;; call resimplify, but this doesn't fix any bugs that I know.
  (setq e (resimplify e))
  ;; Make sure that `e` is in standard representation
  (setq e ($ratdisrep e))

  (residue-by-methods e var p))
 
(defun residue-by-infinity-transform (e x pt)
 "Compute the residue of `e` at infinity by transforming via x = 1/z."
	(and (member pt *infinities*)
	     (residue-by-methods (mul -1 (div (maxima-substitute (div 1 x) x e) (mul x x))) x 0)))

(defun residue-rational (w x pt)
  "Compute the residue of the rational function `w` with respect to variable `x` at point `pt`.
 Arguments:
  `w`  -- a rational expression 
  `x`  -- the variable with respect to which the residue is computed
  `pt` -- the point at which the residue is evaluated

Returns: The residue of the function `x -> w` at `pt`."

  (setq w (sratsimp w))
  (let* ((p ($num w))
         (q ($factor ($denom w)))
         (qfactors (if (mtimesp q) (cdr q) (list q)))
         (cnd)
         (vanish))

    (cond ((ratp w x)
           ;; We assume that the members of qfactors are relatively prime. When the polynomial coefficients are
           ;; rational numbers, this condition is satisfied. But when some coefficients are symbolic, this might
           ;; not be true for all values of the parameters; for these cases, the condition that the members of
           ;; qfactors be relatively prime are

           (let ((cntx ($supcontext)))
             (unwind-protect
                 (progn
                   (dolist (qk qfactors)
                     (dolist (ql qfactors)
                       (assume (ftake '$notequal ($resultant qk ql x) 0))))
                   ;; additionally, assume that the numerator and denominator have no common factors
                   (assume (ftake '$notequal ($resultant p q x) 0))
                   (mtell (intl:gettext "Assuming:  ~M ~%") (fapply 'mand (cdr (mfuncall '$facts))))

                   (catch 'finished
                     (dolist (qk qfactors)
                       (when (eq '$yes ($askequal 0 (square-free-part (maxima-substitute pt x qk))))
                         (setq vanish (square-free-part (maxima-substitute pt x qk)))
                         (let ((n (if (mexptp qk) (third qk) 1)))
                           (throw 'finished
                             (get-taylor-coeff
                               (div p ($first ($divide q (power (sub x pt) n))))
                               x pt (sub n 1) vanish)))))
                     (throw 'finished 0)))
               ($killcontext cntx))))
          (t nil))))

(defun residue-by-taylor (e x pt &optional (n 4) (stop 2))
  "Use a Taylor polynomial to find the residue of `e` at `pt` with respect to `x`. When successful,
  return the residue, otherwise return nil. The fourth argument `n` determines the order of the Taylor
  polynomial. When the order is too small, double `n` and try again. The last argument `stop` counts down
  and ends this recursion when `stop` reaches zero."
  (let ((ee)
        (silent-taylor-flag t)
        ($taylordepth 8)
        ($radexpand nil)
        ($%emode t)
        ($taylor_logexpand nil)
        ($exponentialize nil)
        ($demoivre nil)
        ($taylor_simplifier #'resimplify)
        ($logexpand nil))
    (cond
      ;; This shouldn't happen, but we'll check for this case
      ((member pt *infinities*) nil)

      ((not (tlimp e x)) nil)

      ;; Terminate effort to find the Taylor series when `stop` is zero
      ((eql stop 0) nil)

      (t
       (setq e (logarc-atan2 e))
       (setq ee (catch 'taylor-catch ($taylor e x pt n)))
       (cond
         ;; It is important to check that `taylor` returns a sum of integer powers; here is a case that it
         ;; does not: taylor(log(x),x,0,3) -> log(x)+.... Possibly, we could use taylorinfo to check if
         ;; taylor returns a sum of powers. When taylor was successful, return the residue
         
         ((and ee
               (not (eql ee 0))
               ($polynomialp ee (ftake 'mlist x) #'(lambda (s) (freeof x s)) #'integerp)
               (alike1 ($taylorinfo ee) (ftake 'mlist (ftake 'mlist x pt n))))
          ;; Return the coefficient of 1/(x - pt) from the Taylor expansion
          ($ratdisrep (pscoeff1 ee x -1)))

         ;; Double n and retry if ee was not nil and stop is positive
         ((and ee (> stop 0))
          (residue-by-taylor e x pt (* 2 (max 1 n)) (1- stop)))

         ;; Otherwise, give up
         (t nil))))))

;; I think this code needs a way to detect branch points. It is responsible for the bug
;; residue(1/(sqrt(x^2 - 1)), x, 1)
(defun residue-by-taylor-asym (e x pt &optional (n 4) (stop 2))
  "Use a Taylor polynomial to find the residue of `e` at `pt` with respect to `x`. When successful,
  return the residue, otherwise return nil. The fourth argument `n` determines the order of the Taylor
  polynomial. When the order is too small, double `n` and try again. The last argument `stop` counts down
  and ends this recursion when `stop` reaches zero."
  (let ((ee)
        (silent-taylor-flag t)
        ($taylordepth 8)
        ($radexpand nil)
        ($%emode t)
        ($taylor_logexpand nil)
        ($exponentialize nil)
        ($demoivre nil)
        ($taylor_simplifier #'resimplify)
        ($logexpand nil))
    (cond
      ;; This shouldn't happen, but we'll check for this case
      ((member pt *infinities*) nil)

      ((not (tlimp e x)) nil)

      ;; Terminate effort to find the Taylor series when `stop` is zero
      ((eql stop 0) nil)

      (t
       (setq e (logarc-atan2 e))
       (setq ee (catch 'taylor-catch ($taylor e (ftake 'mlist x pt n '$asymp))))

       (cond
         ;; It is important to check that `taylor` returns a sum of integer powers; here is a case that it
         ;; does not: taylor(log(x),x,0,3) -> log(x)+.... Possibly, we could use taylorinfo to check if
         ;; taylor returns a sum of powers. When taylor was successful, return the residue
         ((and ee
               (not (eql ee 0))
               (alike1 ($taylorinfo ee) (ftake 'mlist (ftake 'mlist x pt n)))
               ($polynomialp ee (ftake 'mlist x) #'(lambda (s) (freeof x s)) #'integerp))

          ;; Return the coefficient of 1/(x - pt) from the Taylor expansion. And yes, for the
          ;; asymp form of taylor expansions, the third argument to `pscoeff1` needs to be 1, 
          ;; not -1.
          (let ((cf ($ratdisrep (pscoeff1 ee x 1))))
            (if (eql cf 0)
                nil
                cf)))

         ;; Double n and retry if ee was not nil and stop is positive
         ((and ee (> stop 0))
          (residue-by-taylor e x pt (* 2 (max 1 n)) (1- stop)))

         ;; Otherwise, give up
         (t nil))))))

(defun residue-nounform (e x pt)
"Construct a symbolic noun form of the residue expression residue(e, x, pt)."
	(ftake ($nounify '$residue) e x pt))

;; This methods doesn't ask questions that some other methods might ask--for example 
;; residue(a/b,x,b).
(defun residue-by-freeof (e x pt)
  "Return 0 if `x` is free of `e`, so the residue is zero. Otherwise, return nil."
  (declare (ignore pt))
  (if (freeof x e)
      0
      nil))

(defun residue-by-simp (e x pt)
 "Uses linearaity to to simplify a residue expression."
  (cond ((mplusp e)
         (fapply 'mplus
                 (mapcar #'(lambda (s) (residue-by-methods s x pt :methods (list 'residue-by-simp))) (cdr e))))
        ((and (mtimesp e) (freeof x (cadr e)))
         (mul (cadr e)
              (residue-by-methods (fapply 'mtimes (cddr e)) x pt :methods (list 'residue-by-simp))))
        (t
         (residue-by-methods e x pt :methods (list 'residue-nounform)))))

(defun residue-by-misc-undefined (e x pt)
"Return a residue nounform when input is nonsence."
  ;; Possibly optionally throw an error instead of a nounform.
  (if (or (mrelationp e) (mrelationp pt) (not ($mapatom x)) (mbagp e) (mbagp pt) (not (freeof x pt)))
     (residue-by-methods e x pt :methods (list 'residue-nounform))
     nil))

(defun residue-by-methods (e x pt &key (methods *residue-methods*))
  "Attempt to compute the residue of expression `e` with respect to variable `x` at point `pt`.

Tries each method in the `methods` list in order. Each method should accept the input `(e x pt)`
and return either a valid residue or nil. Returns the first successful result, or nil if all methods fail.

Optional keyword argument:
  :methods — a list of method symbols to try (default: *residue-methods*)"
  ;; (setq e (mfuncall '$trigsimp e)) ; optional preprocessing
  (let (($gcd '$spmod) ($algebraic t))
  (catch 'finished
    (dolist (fn methods)
      (let ((ans (funcall fn e x pt)))
        (when ans
          (when *residue-method-info*
            (mtell (intl:gettext "Residue method: ~A succeeded.~%") fn))
          (throw 'finished ans)))))))

(defun resm1 (e pt)
   (residue-by-methods e var pt))

(defun resm1-var (x e pt)
	 (residue-by-methods e x pt))

(defun square-free-part (p)
 "Return the square-free part of polynomial `p` using Maxima's $sqfr. Each factor is retained, but 
 its multiplicity is reduced to 1."
  (let ((fp ($sqfr p)))
    (cond
      ((mtimesp fp)
       (fapply 'mtimes (mapcar #'(lambda (s) (if (mexptp s) (second s) s)) (cdr fp))))
      ((mexptp fp)
       (second fp))
      (t fp))))

;; Undone: (a) branch point detector (b) residue by matching

(defun sump (e)
  "Return `t` if `e` is a Maxima sum expression; otherwise return nil."
  (and (consp e) (eq '%sum (caar e))))

;; should I set sumexpand : true & cauchysum : true?  what about
;; residue(exp(1/x)/(x+a),x,0)? I think we need to be more careful with
;; solve. Can solve fail? What about all the solve option variables?
(defun residue-by-powerseries (e x pt)
  "Compute the residue of expression `e` at point `pt` with respect to variable `x`,
   using a power series expansion. Constructs the series in a temporary context
   to avoid assumptions about the summation index, then safely discards it.
   Returns the coefficient of (x - pt)^-1 if identifiable, otherwise NIL."

  ;; When algebraic is true, we get a bit of a mess from powerseries(1/(x^(2/3) + 1),x,0).
  ;; So we'll set algrebraic to false.
  (let ((cntx ($supcontext)) ($sumexpand t) ($cauchysum t) ($algebraic nil) (ps))
    (unwind-protect
        (setq ps (car (errcatch ($intosum ($powerseries e x pt)))))
        ($killcontext cntx))

    (cond
      ((sump ps)
       (let* ((summand (second ps)) ;not correct for iterated sums
              (index (third ps))
              (lo (fourth ps))
              (hi (fifth ps))
              (n (div (mul ($diff summand x) (sub x pt)) summand))
              (nn ($rhs ($first ($solve (ftake 'mequal n -1) index))))
              (cntx ($supcontext)))
         (unwind-protect
         (if (and (eq '$yes ($askinteger nn))
                  (eq '$yes ($ask (ftake 'mleqp lo nn) t))
                  (eq '$yes ($ask (ftake 'mleqp nn hi) t)))
             (coeff (maxima-substitute nn index summand) x -1)
             nil)
         ($killcontext cntx))))
      (t nil))))

;; experimental code--ask a question about a mrelationp expression. When true, assume the fact
;; in the current context.
(defmfun $ask (e &optional (save nil))
  (cond
    ((mrelationp e)
     (let ((ans (ask e)))
       (when save
         (assume e))
       ans))
    (t
     (merror "ask: Expected an <, <=, =, #, >, >= expression, but found ~M ~%" e)
     '$done)))

(defmfun ask (e)
  (let ((answer (mfuncall '$is e)))
    (cond
      ((eq answer t) '$yes)
      ((eq answer nil) '$no)
      (t
       (setq answer (retrieve `((mtext) ,(intl:gettext "Is ") ,e  ,(intl:gettext "? ")) nil))
       (cond
         ((member answer '($no |$n| |$N|) :test #'eq) '$no)
         ((member answer '($yes |$y| |$Y|) :test #'eq) '$yes)
         (t
          (mtell (intl:gettext "Acceptable answers are yes, y, Y, no, n, N. ~%"))
          (ask e)))))))