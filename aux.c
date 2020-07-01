/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics->cn
  *Create Time: 2017-05-30 16:32:17
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aux.h"
#include "str.h"
#include "utils.h"

int type2nbytes[] = {1, 1, 2, 2, 4, 4, 8, 8, 4, 8, 0};

MP_DEF (aux, aux_item_t);

void
aux_item_init (aux_item_t * item)
{
	item->s = str_init ();
}

void
aux_item_free (aux_item_t * item)
{
	str_free (item->s);
}

void
aux_item_copy (aux_item_t * dst, aux_item_t * src)
{
	dst->key[0] = src->key[0];
	dst->key[1] = src->key[1];
	dst->type = src->type;
	dst->deleted = src->deleted;
	dst->tag = src->tag;

	if (src->type == AUX_STR)
		str_copy (dst->s, src->s);
	else
		memcpy (dst->s->s, src->s->s, type2nbytes[src->type]);
}

aux_t *
aux_init (void)
{
	aux_t * a;

	a = (aux_t *) ckalloc (1, sizeof(aux_t));
	a->pool = mp_init (aux, aux_item_init);

	return a;
}

void
aux_clear (aux_t * a)
{
	mp_clear (aux, a->pool, NULL);
}

void
aux_free (aux_t * a)
{
	mp_free (aux, a->pool, aux_item_free);
	free (a);
}

int32_t
aux_cnt (aux_t * a)
{
	return a->pool->n;
}

aux_item_t *
aux_at (aux_t * a, int32_t idx)
{
	return a->pool->pool + idx;
}

aux_item_t *
aux_find (aux_t * a, char * key)
{
	int32_t i;
	aux_item_t * item;

	for (i=0; i<mp_cnt(a->pool); ++i) {
		item = mp_at (aux, a->pool, i);
		if (item->key[0] == key[0]
				&& item->key[1] == key[1]
				&& item->deleted == 0)
			return item;
	}

	return NULL;
}

int
aux_add (aux_t * a, char * key, uint8_t * val, uint8_t type, char tag)
{
	aux_item_t * item;

	item = mp_alloc (aux, a->pool);
	memcpy (item->key, key, 2);
	item->type = type;
	item->deleted = 0;
	item->tag = tag;

	if (type == AUX_STR)
		str_assign (item->s, (const char *)val);
	else
		memcpy (item->s->s, val, type2nbytes[type]);

	return 0;
}

int
aux_del (aux_t * a, char * key)
{
	aux_item_t * item;
	
	if ((item = aux_find(a, key)) == NULL)
		return -1;

	item->deleted = 1;

	return 0;
}

int
aux_edit (aux_t * a, char * key, uint8_t * val)
{
	aux_item_t * item;

	if ((item = aux_find(a, key)) == NULL)
		return -1;

	if (item->type == AUX_STR)
		str_assign (item->s, (const char *)val);
	else
		memcpy (item->s->s, val, type2nbytes[item->type]);

	return 0;
}

int
aux_read (FILE * fp, aux_t * aux)
{
	int32_t i;
	aux_item_t * item;

	fread (&aux->pool->n, 4, 1, fp);
	mp_resize (aux, aux->pool, aux->pool->n);
	for (i=0; i<aux->pool->n; i++) {
		item = mp_at (aux, aux->pool, i);

		fread (item->key, 1, 4, fp);
		fread (&item->tag, 1, 1, fp);
		fread (&item->type, 1, 1, fp);
		fread (&item->deleted, 2, 1, fp);
		if (item->type == AUX_STR) {
			str_read (fp, item->s);
		} else
			fread (item->s->s, 1, 8, fp);
	}

	return 0;
}

int
aux_write (FILE * fp, aux_t * aux)
{
	int32_t i;
	aux_item_t * item;

	fwrite (&aux->pool->n, 4, 1, fp);
	for (i=0; i<mp_cnt(aux->pool); i++) {
		item = mp_at (aux, aux->pool, i);

		fwrite (item->key, 1, 4, fp);
		fwrite (&item->tag, 1, 1, fp);
		fwrite (&item->type, 2, 1, fp);
		fwrite (&item->deleted, 2, 1, fp);
		if (item->type == AUX_STR)
			str_write (fp, item->s);
		else
			fwrite (item->s->s, 8, 1, fp);
	}

	return 0;
}

int
aux_copy (aux_t * dst, aux_t * src)
{
	mp_copy (aux, dst->pool, src->pool, aux_item_copy);

	return 0;
}

int
aux2hts_aux (aux_item_t * item, str_t * s)
{
	int32_t l_item;

	if (item->type == AUX_STR) {
		l_item = 2 + 1 + item->s->l+1;
		str_resize (s, s->l+l_item);
		memcpy (s->s+s->l, item->key, 2); s->l += 2;
		s->s[s->l] = item->tag; s->l += 1;
		memcpy (s->s+s->l, item->s->s, item->s->l+1);
		s->l += item->s->l+1;
	} else {
		l_item = 2 + 1 + type2nbytes[item->type];
		str_resize (s, s->l+l_item);
		memcpy (s->s+s->l, item->key, 2); s->l += 2;
		s->s[s->l] = item->tag; s->l += 1;
		memcpy (s->s+s->l, item->s->s, type2nbytes[item->type]);
		s->l += type2nbytes[item->type];
	}

	return 0;
}
