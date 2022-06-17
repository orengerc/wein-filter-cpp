#include "structs.h"

/******************************************/
/*        FUNCTIONS DECLARATIONS          */
/******************************************/

void free_vec(Vec **p_vec);

/***********************************************/
/*        H FUNCTIONS IMPLEMENTATIONS          */
/***********************************************/

Timeline *alloc_time_line()
{
  Timeline *timeline = calloc (1, sizeof(Timeline));
  return timeline;
}

TimeState *alloc_time_state(long double t,Vec *r, Vec *v, Vec *a)
{
  TimeState *time_state = malloc (sizeof (TimeState));
  if (!time_state)
  {
    return NULL;
  }
  time_state->time = t;
  time_state->r = r;
  time_state->v = v;
  time_state->a = a;
  time_state->next = NULL;
  return time_state;
}

Vec *alloc_vec(long double y, long double z)
{
  Vec *vec = malloc (sizeof (Vec));
  if (!vec) {return NULL;}
  vec->_y = y;
  vec->_z = z;
  return vec;
}

void free_time_line(Timeline **p_timeline)
{
  Timeline *timeline = *p_timeline;
  if (!timeline) {return;}
  TimeState *next_time_state = timeline->first;
  for (int i = 0; i < timeline->size; ++i)
  {
    TimeState *curr_time_state = next_time_state;
    next_time_state = curr_time_state->next;
    free_time_state(&curr_time_state);
  }
  free(timeline);
  *p_timeline = NULL;
}

void mult_vec_scalar(Vec *vec, long double n)
{
  if (!vec) {return;}
  vec->_y = vec->_y * n;
  vec->_z = vec->_z * n;
}

void add_to_vec_scalar(Vec *vec, long double n)
{
  if (!vec) {return;}
  vec->_y = vec->_y + n;
  vec->_z = vec->_z + n;
}

void mult_vec_vec(Vec *first, Vec *second)
{
  if (!first || !second){return;}
  first->_y = first->_y * second->_y;
  first->_z = first->_z * second->_z;
}

void add_to_vec_vec(Vec *first, Vec *second)
{
  if (!first || !second){return;}
  first->_y = first->_y + second->_y;
  first->_z = first->_z + second->_z;
}

TimeState *clone_time_state(TimeState *time_state)
{
  Vec *r = alloc_vec (time_state->r->_y, time_state->r->_z);
  Vec *v = alloc_vec (time_state->v->_y, time_state->v->_z);
  Vec *a = alloc_vec (time_state->a->_y, time_state->a->_z);

  if (!r || !v || !a)
  {
    free (r);
    free (v);
    free (a);
    return NULL;
  }

  return alloc_time_state (time_state->time, r, v, a);
}

long double abs_d(long double val)
{
  if (val < 0)
  {
    return (-1) * val;
  }
  return val;
}

/***************************/
/*        HELPERS          */
/***************************/

void free_time_state(TimeState **p_time_state)
{
  TimeState *time_state = *p_time_state;
  if (!time_state) {return;}
  free_vec(&time_state->r);
  free_vec(&time_state->v);
  free_vec(&time_state->a);
  free (time_state);
  *p_time_state = NULL;
}

void free_vec(Vec **p_vec)
{
  free(*p_vec);
  *p_vec = NULL;
}

