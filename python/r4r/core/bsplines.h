
static PyObject* cox_de_boor(PyObject* self, PyObject* args);

static PyObject* find_knot_span(PyObject* self, PyObject* args);

static PyObject* evaluate_curve(PyObject* self, PyObject* args);

static PyObject* evaluate_surface(PyObject* self, PyObject* args);

static PyObject* unwrap_phase(PyObject* self, PyObject* args);

static PyObject* characteristic_function(PyObject* self, PyObject* args);

/*!
 * \brief Cox-de-Boor recursion.
 * \param[in] order polynomial degree of the spline plus 1
 * \param[in] knot knot vector
 * \param[in] t location at which to evaluate the B-spline basis
 * \param[in] N values of the basis at \f$t\f$
 */
void _cox_de_boor(int order, double* knot, double t, double* N);

/*!
 * \brief Cox-de-Boor recursion including derivatives.
 * \param[in] order polynomial degree of the spline plus 1
 * \param[in] knot knot vector
 * \param[in] t location at which to evaluate the B-spline basis ans its derivatives
 * \param[in] N values of the basis and its derivaties at \f$t\f$
 */
void _cox_de_boor_derivatives(int order, double* knot, int der_count, double* N);

/*!
 * \brief Binary search of knot vector for knot span.
 * \param[in] t key
 * \param[in] knots knot vector
 * \param[in] lower lower bound on range
 * \param[in] upper upper bound on range
 * \return knot span
 *
 * Make sure before calling the function that t is not the last knot. The function
 * should then be called with upper=number of knots - 2, as the items of the recursions
 * are the intervals.
 *
 */
int _find_knot_span(double t, double* knots, int lower, int upper);

void _unwrap_phase(double* x, size_t n);
