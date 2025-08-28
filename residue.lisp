(in-package :maxima)

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
 of the solutions and a CL list of the corresponding multiplicities. When solve fails, 
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
(defmfun $residue (e x pt)
  ;; Error when var isn't a mapatom
  (when (not ($mapatom x))
    (merror (intl:gettext "residue: second argument must be a mapatom; found ~M") x))
  ;; Error when `p` depends on `var`
  (when (not (freeof x pt))
    (merror (intl:gettext "residue: third argument must not depend on second argument; found ~M") pt))
  ;; Error when `e` is an mbag
  (when (mbagp e)
	(merror (intl:gettext "residue: first argument must not be an mbag; found ~M") e))
   ;; Error when `pt` is an mbag
  (when (mbagp pt)
	(merror (intl:gettext "residue: third argument must not be an mbag; found ~M") pt))

  ;; Unfortunately, when running the testsuite, sometimes `e` is not simplified. So, let us
  ;; call resimplify, but this doesn't fix any bugs that I know.
  (setq e (resimplify e))
  ;; Make sure that `e` is in standard representation
  (setq e ($ratdisrep e))

  (residue-by-methods e x pt))
 
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

  ;; When the denominator of `e` is a polynomial, ask the user if each factor vanishes at `pt`.  When there is
  ;; a vanishing factor, set taylor_simplifier to ratsubst zero for this factor.
  (let ((Q ($denom e)) (vanish nil) (cnd))
    (when ($polynomialp Q (ftake 'mlist x) #'(lambda (s) (freeof x s)))
      (setq Q (square-free-part Q))
      (setq Q (if (mtimesp Q) (cdr Q) (list Q)))
      (catch 'finished
        (dolist (qk Q)
          (setq cnd ($askequal (maxima-substitute pt x qk) 0))
          (when (eq cnd '$yes)
            (setq vanish (maxima-substitute pt x qk))
            (throw 'finished t)))))

    (let ((ee)
          (silent-taylor-flag t)
          ($taylordepth 8)
          ($radexpand nil)
          ($%emode t)
          ($taylor_logexpand nil)
          ($exponentialize nil)
          ($demoivre nil)
          ($taylor_simplifier (if vanish
                                  #'(lambda (s) ($ratsubst 0 vanish s))
                                  #'resimplify))
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
           (t nil)))))))

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
        ;($exponentialize nil)
        ;($demoivre nil)
        ($errormsg nil)
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
       ;; errcatch works around the taylor(exp(x)/(x-a)^n,[x,a,4,asymp]) -> Lisp error bug
       (setq ee (car (errcatch (catch 'taylor-catch ($taylor e (ftake 'mlist x pt n '$asymp))))))
       (cond
         ;; It is important to check that `taylor` returns a sum of integer powers; here is a case that it
         ;; does not: taylor(log(x),x,0,3) -> log(x)+...I'm not sure the taylorinfo check is needed, but
         ;; note that in Maxima 5.48, taylorinfo ignores the asymp identitfer, so using taylorinfo this
         ;; way requires a fix to hayat.
         ((and ee
               (not (eql ee 0))
               (alike1 ($taylorinfo ee) (ftake 'mlist (ftake 'mlist x pt n '$asymp)))
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
  (list (list '$residue 'simp) e x pt))

;; We include this method because it doesn't ask questions that some other methods might ask--for example 
;; for residue(a/b,x,b) some other methods might ask if b is zero.
(defun residue-by-freeof (e x pt)
  "Return 0 if `x` is free of `e`; otherwise, return nil."
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
"Return a residue nounform when input is an inequation, mbag, the variable isn't a mapatom, or the 
 residue point depends on the variable; otherwise, return nil."
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
  ;(setq e (mfuncall '$trigsimp e)) ; optional preprocessing

  ;; Using sqrtdenest on pt fixes some bad errors including integrate((-14*x^2-32)/(x^4+3*x^2+1)^2,x,0,inf);
  (setq pt ($sqrtdenest pt))
  (let (($gcd '$spmod) ($algebraic t))
  (catch 'finished
    (dolist (fn methods)
      (let ((ans (funcall fn e x pt)))
        (when ans
          (when *residue-method-info*
            (mtell (intl:gettext "Residue method: ~A succeeded.~%") fn))
          (throw 'finished ans)))))))

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
  ;; So we'll locally set algebraic to false.
  (let ((ps nil))
    (let ((cntx ($supcontext)) ($sumexpand nil) ($cauchysum nil) ($algebraic nil))
      (unwind-protect
          (setq ps ($intosum (car (errcatch ($powerseries e x pt)))))
        ($killcontext cntx)))
    (cond
      ;; when the powerseries involves an %at expression, the series is likely 
      ;; wrong, return nil.
      ((not (freeof '%at ps)) nil)
      ((and (sump ps) (freeof '%sum ($args ps)))
       (let* ((summand (second ps)) ;incorrect for iterated sums
              (index (third ps))
              (lo (fourth ps))
              (hi (fifth ps))
              (n (div (mul ($diff summand x) (sub x pt)) summand))
              (nn ($rhs ($first ($solve (ftake 'mequal n -1) index))))
              (cntx ($supcontext)))
         (unwind-protect
             (if (and (eq '$yes ($askinteger nn))
                      (eq '$yes ($ask_relational (ftake 'mleqp lo nn) t))
                      (eq '$yes ($ask_relational (ftake 'mleqp nn hi) t)))
                 (coeff (maxima-substitute nn index summand) x -1)
               nil)
           ($killcontext cntx))))
      (t nil))))

;; experimental code--ask a question about a mrelationp expression. When true, assume the fact
;; in the current context.  Possibly "ask" implies that this function does more than it does--it
;; doesn't allow, for example ask(integerp(zzz)). Maybe it either needs to be extended, or the 
;; name needs to be changed. 
(defmfun $ask_relational (e &optional (remember t))
  (cond
    ((mrelationp e)
     (let ((ans (ask-relational-helper e)))
       (when (eq t remember)
         (assume e))
       ans))
    (t
     (merror (intl:gettext "ask: Expected a relational expresion (<, <=, =, #, >, >=), but got ~M ~%") e)
     nil)))

;; Bugs & things to thing about:  
;; (a) ask(x # 3) is accepted, but assume(x # 3) is not valid. I'm not sure what I want.
;; (c) I'm not sure that (mfuncall '$is e) is what I want?
(defmfun ask-relational-helper (e)
  (let ((answer (mfuncall '$is e)))
    (cond
      ((eq answer t) '$yes)
      ((eq answer nil) '$no)
      (t
       (setq answer (string-downcase (symbol-name 
               (retrieve `((mtext) ,(intl:gettext "Is ") ,e  ,(intl:gettext "? ")) nil))))
       (cond
         ((member answer '("$yes" "$y") :test #'string=) '$yes)
         ((member answer '("$no" "$n") :test #'string=) '$no)
         (t
          (mtell (intl:gettext "Acceptable answers are yes, y, no, n (case independent). ~%"))
          (ask-relational-helper e)))))))

;; with CCL, but *not* SBCL, running rtest_residue gives the error: Odd-length property list in REMF.
;; Wrapping the call to remf in zl-remprop with ignore-errors allows the tests to run to completion.
;; This code tries to find out what is going on.
(defvar *yikes* nil)
(defun zl-remprop (sym indicator)
  (if (symbolp sym)
      (remprop sym indicator)
    (unless (atom sym)
      (multiple-value-bind (result condition)
          (ignore-errors
            (remf (cdr sym) indicator))
        (when condition
		   (mtell "remf error:  result = ~M ; condition = ~M ~%" result condition)
		   (let ((*print-circle* t)) (print (cdr sym)))
		   (let ((*print-circle* t)) (print `(sym = ,sym)))
		   (let ((*print-circle* t)) (print `(IND = ,indicator)))
          (push (ftake 'mlist result condition (list-length (cdr sym))) *yikes*))))))

;; Unless you build Maxima from the curent source, there is a bug in taylor info. Here is 
;; a fix for that bug.
(defun taylor-info (q)
  (let ((acc-var nil) (acc-pt nil) (acc-ord nil) (qk) (acc))
    (cond ((null q) nil)
	  (t
	   (setq qk (pop q))
	   (cond ((and (fourth qk) (consp (fourth qk)) (eq (caar (fourth qk)) 'multivar)) nil)
		 ((and (fourth qk) (consp (fourth qk)) (eq (caar (fourth qk)) 'multi))
		  (while (and (fourth qk) (consp (fourth qk)) (eq (caar (fourth qk)) 'multi))
		    (setq acc nil)
		    (push (taylor-trunc qk) acc-ord)
		    (push (exp-pt qk) acc-pt)
		    (push (datum-var qk) acc-var)
		    (setq qk (pop q)))
		  (push '(mlist) acc-ord)
		  (push '(mlist) acc-pt)
		  (push '(mlist) acc-var)
		  (setq q (taylor-info q))
		  (if (null q) (list acc-var acc-pt acc-ord) (append q (list acc-var acc-pt acc-ord))))

		 (t
		  (setq acc (if (and (fourth qk) (consp (fourth qk)) (eq '$asymp (caar (fourth qk))))
				(list '$asymp) nil))
		  (push (taylor-trunc qk) acc)
		  (push (exp-pt qk) acc)
		  (push (datum-var qk) acc)
		  (push '(mlist) acc)
		  (setq q (taylor-info q))
		  (if (null q) (list acc) (append q (list acc)))))))))