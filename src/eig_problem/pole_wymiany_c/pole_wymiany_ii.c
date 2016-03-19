#include <stdio.h>
#include <complex.h>
#include "uthash.h"

typedef struct vec2d {
	long long x_;
	long long y_;
} vec2d;

typedef struct vec2d_float {
	double x_;
	double y_;
} vec2d_float;

/* struct used as key for
 * dlugosc_wymiany and magnetyzacja values
 * */
typedef struct vec_key {
	long long x_;
	long long y_;
} vec_key;

/* struct for caching dlugosc_wymiany values */
typedef struct dw {
	vec_key key_;
	double complex z_;
	UT_hash_handle hh;
} dw;

/* struct for caching magnetyzacja values */
typedef struct m {
	vec_key key_;
	double complex z_;
	UT_hash_handle hh;
} m;

dw* dlugosc_wymiany = NULL;
m* magnetyzacja = NULL;

int liczba_wektorow = 0;
vec2d* lista_wektorow = NULL;
double H0 = 0.0;
double size = 0.0;

void add_dlugosc_wymiany_value(long long x, long long y, double rz, double iz)
{
	dw* dwv;
	vec_key dwk = {x, y};
	HASH_FIND(hh, dlugosc_wymiany, &dwk, sizeof(vec_key), dwv);
	if (dwv == NULL) {
		dw* new_dw = malloc(sizeof(dw));
		new_dw->key_.x_ = x;
		new_dw->key_.y_ = y;
		new_dw->z_ = rz + iz;
#ifdef DEBUG
		fprintf(stderr, "[c] dlugosc_wymiany (%lld, %lld) = (%e, %e)\n", x, y, rz, iz);
#endif
		HASH_ADD(hh, dlugosc_wymiany, key_, sizeof(vec_key), new_dw);
	}
}

void add_magnetyzacja_value(long long x, long long y, double rz, double iz)
{
	m* mv;
	vec_key mk = {x, y};
	HASH_FIND(hh, magnetyzacja, &mk, sizeof(vec_key), mv);
	if (mv == NULL) {
		m* new_m = malloc(sizeof(m));
		new_m->key_.x_ = x;
		new_m->key_.y_ = y;
		new_m->z_ = rz + iz;
#ifdef DEBUG
		fprintf(stderr, "[c] magnetyzacja (%lld, %lld) = (%e, %e)\n", x, y, rz, iz);
#endif
		HASH_ADD(hh, magnetyzacja, key_, sizeof(vec_key), new_m);
	}
}

void init_lista_wektorow(vec2d* lw, int n, double h0, double a)
{
	liczba_wektorow = n;
	lista_wektorow = malloc(sizeof(vec2d) * n);
	H0 = h0;
    size = a;
	for (int i = 0; i < n; ++i) {
		lista_wektorow[i].x_ = lw[i].x_;
		lista_wektorow[i].y_ = lw[i].y_;
#ifdef DEBUG
		fprintf(stderr, "[c] initializing lista_wektorow[%d] = (%lld, %lld)\n",
			i, lw[i].x_, lw[i].y_);
		fprintf(stderr, "[c] set_value lista_wektorow[%d] = (%lld, %lld)\n",
			i,
			lista_wektorow[i].x_,
			lista_wektorow[i].y_);
#endif
	}
}

void tmp_value(vec2d* w1, vec2d* w2, vec2d_float* wq, vec2d_float* result_tmp)
{
	double complex tmp = 0.0 + 0.0;
#ifdef BTMPV
	fprintf(stderr, "[c] w1 = (%lld, %lld)\n", w1->x_, w1->y_);
	fprintf(stderr, "[c] w2 = (%lld, %lld)\n", w2->x_, w2->y_);
	fprintf(stderr, "[c] wq = (%f, %f)\n", wq->x_, wq->y_);
#endif
	for (int i = 0; i < liczba_wektorow; ++i) {
		dw* dwv;
		m* mv;
		vec_key dwk = {lista_wektorow[i].x_ - w2->x_,
			       lista_wektorow[i].y_ - w2->y_};
		vec_key mk = {w1->x_ - lista_wektorow[i].x_,
			      w1->y_ - lista_wektorow[i].y_};
		HASH_FIND(hh, dlugosc_wymiany, &dwk, sizeof(vec_key), dwv);
		HASH_FIND(hh, magnetyzacja, &mk, sizeof(vec_key), mv);
#ifdef BTMPV
		fprintf(stderr, "[c] H0 = %f\n", H0);
		fprintf(stderr, "[c] vec_l = (%lld, %lld)\n",
				lista_wektorow[i].x_,
				lista_wektorow[i].y_);
		fprintf(stderr, "[c] sdw (%lld, %lld) = (%e, %e)\n",
			lista_wektorow[i].x_ - w2->x_,
			lista_wektorow[i].y_ - w2->y_,
			creal(dwv->z_),
			cimag(dwv->z_)
			);
		fprintf(stderr, "[c] sm (%lld, %lld) = (%e, %e)\n",
			w1->x_ - lista_wektorow[i].x_,
			w1->y_ - lista_wektorow[i].y_,
			creal(mv->z_),
			cimag(mv->z_)
			);
#endif
		tmp += ((wq->x_ + w2->x_ / size) * (wq->x_ + lista_wektorow[i].x_ / size)
		      + (wq->y_ + w2->y_ / size) * (wq->y_ + lista_wektorow[i].y_ / size))
			* dwv->z_ * mv->z_ / H0;
#ifdef TMPV
		fprintf(stderr, "[c] inside tmp = (%f, %f)\n", creal(tmp), cimag(tmp));
#endif
	}
#ifdef TMPV
	fprintf(stderr, "[c] outside tmp = (%f, %f)\n", creal(tmp), cimag(tmp));
#endif
	result_tmp->x_ = creal(tmp);
	result_tmp->y_ = cimag(tmp);
}
