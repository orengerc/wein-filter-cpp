#include "methods.h"

/******************************************/
/*        FUNCTIONS DECLARATIONS          */
/******************************************/

void free_and_NULLify (Vec **vec);
Vec *get_analytic_r (long double t);
Vec *get_analytic_v (long double t);
Vec *get_analytic_a (Vec *v);
Vec *euler_method_helper (Vec *val, Vec *der, long double Dt);
Vec *get_a (Vec *v);
Vec *get_k_1_v (Vec *_v, long double Dt);
Vec *get_k_2_v (Vec *_v, Vec *k_1, long double Dt);
Vec *get_k_1_r (Vec *_v, long double Dt);
Vec *get_k_2_r (Vec *_v, Vec *k_1, long double Dt);
Vec *get_midpoint_vec (Vec *vec, Vec *k_2);
void free_midpoint (Vec **a, Vec **v, Vec **r, Vec **k_1_v, Vec **k_2_v, Vec
**k_1_r, Vec **k_2_r);
Vec *get_k_3_v (Vec *_v, Vec *k_2, long double Dt);
Vec *get_k_4_v (Vec *_v, Vec *k_3, long double Dt);
Vec *get_k_3_r (Vec *_v, Vec *k_2, long double Dt);
Vec *get_k_4_r (Vec *_v, Vec *k_3, long double Dt);
Vec *get_runge_kutta_vec (Vec *vec, Vec *k_1, Vec *k_2, Vec *k_3, Vec *k_4);
void free_runge_kutta (Vec **a, Vec **v, Vec **r, Vec **k_1_v, Vec **k_2_v,
                       Vec **k_3_v, Vec **k_4_v, Vec **k_1_r, Vec **k_2_r,
                       Vec **k_3_r, Vec **k_4_r);

/***********************************************/
/*        H FUNCTIONS IMPLEMENTATIONS          */
/***********************************************/

TimeState *analytic_method (TimeState *curr_time_state, long double Dt)
{
  long double _time = curr_time_state->time;

  Vec *r = get_analytic_r (_time + Dt);
  Vec *v = get_analytic_v (_time + Dt);
  Vec *a = get_analytic_a (v);

  TimeState *next_time_state = alloc_time_state (_time + Dt,
                                                 r, v, a);
  if (!next_time_state || !a || !v || !r)
  {
    free (a);
    free (v);
    free (r);
    free (next_time_state);
    return NULL;
  }
  return next_time_state;
}

TimeState *euler_method (TimeState *curr_time_state, long double Dt)
{
  Vec *_r = curr_time_state->r;
  Vec *_v = curr_time_state->v;
  Vec *_a = curr_time_state->a;

  Vec *a = get_a (_v);
  Vec *v = euler_method_helper (_v, _a, Dt);
  Vec *r = euler_method_helper (_r, _v, Dt);

  TimeState *next_time_state = alloc_time_state (curr_time_state->time + Dt,
                                                 r, v, a);
  if (!next_time_state || !a || !v || !r)
  {
    free (a);
    free (v);
    free (r);
    free (next_time_state);
    return NULL;
  }
  return next_time_state;
}

TimeState *midpoint_method (TimeState *curr_time_state, long double Dt)
{
  long double _time = curr_time_state->time;
  Vec *_v = curr_time_state->v;
  Vec *_r = curr_time_state->r;

  Vec *a = get_a (_v);

  Vec *k_1_v = get_k_1_v (_v, Dt);
  Vec *k_2_v = get_k_2_v (_v, k_1_v, Dt);
  Vec *v = get_midpoint_vec (_v, k_2_v);

  Vec *k_1_r = get_k_1_r (_v, Dt);
  Vec *k_2_r = get_k_2_r (_v, k_1_v, Dt);
  Vec *r = get_midpoint_vec (_r, k_2_r);

  if (!r || !v || !a)
  {
    free_midpoint (&a, &v, &r, &k_1_v, &k_2_v, &k_1_r, &k_2_r);
    return NULL;
  }
  TimeState *next_time_state = alloc_time_state (_time + Dt, r, v, a);
  free_midpoint (NULL, NULL, NULL, &k_1_v, &k_2_v, &k_1_r, &k_2_r);
  return next_time_state;
}

TimeState *runge_kutta_method (TimeState *curr_time_state, long double Dt)
{
  long double _time = curr_time_state->time;
  Vec *_r = curr_time_state->r;
  Vec *_v = curr_time_state->v;

  Vec *a = get_a (_v);

  Vec *k_1_v = get_k_1_v (_v, Dt);
  Vec *k_2_v = get_k_2_v (_v, k_1_v, Dt);
  Vec *k_3_v = get_k_3_v (_v, k_2_v, Dt);
  Vec *k_4_v = get_k_4_v (_v, k_3_v, Dt);
  Vec *v = get_runge_kutta_vec (_v, k_1_v, k_2_v, k_3_v, k_4_v);

  Vec *k_1_r = get_k_1_r (_v, Dt);
  Vec *k_2_r = get_k_2_r (_v, k_1_v, Dt);
  Vec *k_3_r = get_k_3_r (_v, k_2_v, Dt);
  Vec *k_4_r = get_k_4_r (_v, k_3_v, Dt);
  Vec *r = get_runge_kutta_vec (_r, k_1_r, k_2_r, k_3_r, k_4_r);

  if (!a || !v || !r)
  {
    free_runge_kutta (&a, &v, &r, &k_1_v, &k_2_v, &k_3_v, &k_4_v, &k_1_r,
                      &k_2_r, &k_3_r, &k_4_r);
    return NULL;
  }
  TimeState *next_time_state = alloc_time_state (_time + Dt, r, v, a);
  free_runge_kutta (NULL, NULL, NULL, &k_1_v, &k_2_v, &k_3_v, &k_4_v, &k_1_r,
                    &k_2_r, &k_3_r, &k_4_r);
  return next_time_state;
}

NEXT_STEP_METHOD *get_method (Method method)
{
  switch (method)
  {
    case ANALYTIC:
      return analytic_method;

    case EULER:
      return euler_method;

    case MIDPOINT:
      return midpoint_method;

    case RUNGE_KUTTA:
      return runge_kutta_method;
  }
  return NULL;
}

/***********************************/
/*        GENERAL HELPERS          */
/***********************************/

void free_and_NULLify (Vec **vec)
{
  if (!vec)
  {
    return;
  }
  free (*vec);
  *vec = NULL;
}

Vec *get_a (Vec *v)
{
  if (!v)
  { return NULL; }
  long double w = (q / m) * B;
  return alloc_vec (((q * E) / m) - (w * v->_z), w * v->_y);
}

/************************************/
/*        ANALYTIC HELPERS          */
/************************************/

Vec *get_analytic_r (long double t)
{
  long double w = (q / m) * B;
  long double y = ((2 * E) / (w * B)) * (cosl (w * t) - 1);
  long double z = ((2 * E) / (w * B)) * sinl (w * t) + (E / B) * t;
  return alloc_vec (y, z);
}

Vec *get_analytic_v (long double t)
{
  long double w = (q / m) * B;
  long double y = (E / B) * (-2 * sinl (w * t));
  long double z = (E / B) * (2 * cosl (w * t) + 1);
  return alloc_vec (y, z);
}

Vec *get_analytic_a (Vec *v)
{
  long double y = (q / m) * (E - B * v->_z);
  long double z = (q / m) * (B * v->_y);
  return alloc_vec (y, z);
}

/*********************************/
/*        EULER HELPERS          */
/*********************************/

Vec *euler_method_helper (Vec *val, Vec *der, long double Dt)
{
  double y = val->_y + der->_y * Dt;
  double z = val->_z + der->_z * Dt;
  return alloc_vec (y, z);
}

/************************************/
/*        MIDPOINT HELPERS          */
/************************************/

Vec *get_k_1_v (Vec *_v, long double Dt)
{
  if (!_v)
  { return NULL; }
  Vec *k_1 = get_a (_v);
  mult_vec_scalar (k_1, Dt);
  return k_1;
}

Vec *get_k_2_v (Vec *_v, Vec *k_1, long double Dt)
{
  if (!_v || !k_1)
  { return NULL; }
  Vec *tmp = alloc_vec (k_1->_y, k_1->_z);
  mult_vec_scalar (tmp, 0.5);
  add_to_vec_vec (tmp, _v);
  Vec *k_2 = get_a (tmp);
  if (!k_2)
  { return NULL; }
  mult_vec_scalar (k_2, Dt);
  free (tmp);
  return k_2;
}

Vec *get_k_1_r (Vec *_v, long double Dt)
{
  if (!_v)
  { return NULL; }
  long double r_y = _v->_y * Dt;
  long double r_z = _v->_z * Dt;
  return alloc_vec (r_y, r_z);
}

Vec *get_k_2_r (Vec *_v, Vec *k_1, long double Dt)
{
  if (!_v || !k_1)
  { return NULL; }
  long double r_y = (_v->_y + (k_1->_y * 0.5)) * Dt;
  long double r_z = (_v->_z + (k_1->_z * 0.5)) * Dt;
  return alloc_vec (r_y, r_z);
}

Vec *get_midpoint_vec (Vec *vec, Vec *k_2)
{
  if (!vec || !k_2)
  { return NULL; }
  long double y = vec->_y + k_2->_y;
  long double z = vec->_z + k_2->_z;
  return alloc_vec (y, z);
}

void free_midpoint (Vec **a, Vec **v, Vec **r, Vec **k_1_v, Vec **k_2_v, Vec
**k_1_r, Vec **k_2_r)
{
  free_and_NULLify (a);
  free_and_NULLify (v);
  free_and_NULLify (r);
  free_and_NULLify (k_1_v);
  free_and_NULLify (k_2_v);
  free_and_NULLify (k_1_r);
  free_and_NULLify (k_2_r);
}

/***************************************/
/*        RUNGE KUTTA HELPERS          */
/***************************************/

Vec *get_k_3_v (Vec *_v, Vec *k_2, long double Dt)
{
  return get_k_2_v (_v, k_2, Dt);
}

Vec *get_k_4_v (Vec *_v, Vec *k_3, long double Dt)
{
  Vec *tmp = alloc_vec (k_3->_y, k_3->_z);
  add_to_vec_vec (tmp, _v);
  Vec *k_4 = get_a (tmp);
  if (!k_4)
  { return NULL; }
  mult_vec_scalar (k_4, Dt);
  free (tmp);
  return k_4;
}

Vec *get_k_3_r (Vec *_v, Vec *k_2, long double Dt)
{
  return get_k_2_r (_v, k_2, Dt);
}

Vec *get_k_4_r (Vec *_v, Vec *k_3, long double Dt)
{
  if (!_v || !k_3)
  { return NULL; }
  long double r_y = (_v->_y + k_3->_y) * Dt;
  long double r_z = (_v->_z + k_3->_z) * Dt;
  return alloc_vec (r_y, r_z);
}

Vec *get_runge_kutta_vec (Vec *vec, Vec *k_1, Vec *k_2, Vec *k_3, Vec *k_4)
{
  if (!vec || !k_1 || !k_2 || !k_3 || !k_4)
  { return NULL; }
  long double y = vec->_y + (1.f / 6) * (k_1->_y + 2 * k_2->_y + 2 * k_3->_y +
                                    k_4->_y);
  long double z = vec->_z + (1.f / 6) * (k_1->_z + 2 * k_2->_z + 2 * k_3->_z +
                                    k_4->_z);
  return alloc_vec (y, z);
}

void free_runge_kutta (Vec **a, Vec **v, Vec **r, Vec **k_1_v, Vec **k_2_v,
                       Vec **k_3_v, Vec **k_4_v, Vec **k_1_r, Vec **k_2_r,
                       Vec **k_3_r, Vec **k_4_r)
{
  free_and_NULLify (a);
  free_and_NULLify (v);
  free_and_NULLify (r);
  free_and_NULLify (k_1_v);
  free_and_NULLify (k_2_v);
  free_and_NULLify (k_3_v);
  free_and_NULLify (k_4_v);
  free_and_NULLify (k_1_r);
  free_and_NULLify (k_2_r);
  free_and_NULLify (k_3_r);
  free_and_NULLify (k_4_r);
}