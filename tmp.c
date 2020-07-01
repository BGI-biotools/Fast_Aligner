/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-08-29 15:32:44
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "str.h"
#include "tmp.h"
#include "mp.h"
#include "utils.h"

#define FILE_MAX 64
#define BIT4FILE_MAX 6
#define SDIR_MAX 64

typedef struct {
	str_t * path;
	str_set_t * files;
	int32_t file_cnt;
	int32_t sdir_cnt;
} sub_dir_t;

MP_DEF (sub_dir, sub_dir_t);

struct tmp_dir_s {
	str_t * root_path;
	mp_t(sub_dir) * sub_dirs;
	pthread_mutex_t mtx;
	sub_dir_t * cur_sub_dir;
	int32_t sub_dir_idx;
};

static void
sub_dir_init (sub_dir_t * sub_dir)
{
	sub_dir->path = str_init ();
	sub_dir->files = str_set_init ();
	sub_dir->file_cnt = 0;
	sub_dir->sdir_cnt = 0;
}

static void
sub_dir_free (sub_dir_t * sub_dir)
{
	int32_t i;
	str_t * file;

	for (i=0; i<str_set_cnt(sub_dir->files); ++i) {
		file = str_set_at (sub_dir->files, i);
	}
	str_set_free (sub_dir->files);
	str_free (sub_dir->path);
}

static str_t *
next_valid_tmp_file (tmp_dir_t * tmp_dir, str_t ** dir, int32_t * file_idx)
{
	sub_dir_t * sdir;
	sub_dir_t * pdir;

	sdir = tmp_dir->cur_sub_dir;
	if (sdir->file_cnt < FILE_MAX) {
		*file_idx = sdir->file_cnt++;
		*dir = sdir->path;
		return str_set_alloc (sdir->files);
	}

	++tmp_dir->sub_dir_idx;
	sdir = mp_alloc (sub_dir, tmp_dir->sub_dirs);
	pdir = mp_at (sub_dir, tmp_dir->sub_dirs, (tmp_dir->sub_dir_idx-1)>>BIT4FILE_MAX);
	str_resize (sdir->path, pdir->path->l+32);
	sprintf (sdir->path->s, "%s/sdir_%d", pdir->path->s, pdir->sdir_cnt++);
  sdir->path->l = strlen (sdir->path->s);
	ckcreate_dir (sdir->path->s);
	tmp_dir->cur_sub_dir = sdir;

	*file_idx = sdir->file_cnt++;
	*dir = sdir->path;
	return str_set_alloc (sdir->files);
}

tmp_dir_t *
tmp_dir_init (const char * dir)
{
	int l_str;
	tmp_dir_t * tmp_dir;

	tmp_dir = (tmp_dir_t *) ckalloc (1, sizeof(tmp_dir_t));
	tmp_dir->root_path = str_init ();
  str_assign (tmp_dir->root_path, dir);
	dir = tmp_dir->root_path->s;
	ckcreate_dir (dir);

	tmp_dir->sub_dirs = mp_init (sub_dir, sub_dir_init);
	tmp_dir->cur_sub_dir = mp_alloc (sub_dir, tmp_dir->sub_dirs);
	tmp_dir->sub_dir_idx = 0;
	str_assign (tmp_dir->cur_sub_dir->path, dir);

	pthread_mutex_init (&tmp_dir->mtx, NULL);

	return tmp_dir;
}

void
tmp_dir_free (tmp_dir_t * tmp_dir)
{
	mp_free (sub_dir, tmp_dir->sub_dirs, sub_dir_free);
	//rm1dir (tmp_dir->root_path->s);
	str_free (tmp_dir->root_path);
	free (tmp_dir);
}

char *
tmp_dir_add1file (tmp_dir_t * tmp_dir, const char * surfix)
{
	str_t * s;
	str_t * dir;
	int32_t file_idx;

	pthread_mutex_lock (&tmp_dir->mtx);
	s = next_valid_tmp_file (tmp_dir, &dir, &file_idx);
	pthread_mutex_unlock (&tmp_dir->mtx);

  if (surfix == NULL) {
	  str_resize (s, dir->l+64);
	  sprintf (s->s, "%s/%d.tmpf", dir->s, file_idx);
  } else {
	  str_resize (s, dir->l+64+strlen(surfix));
    sprintf (s->s, "%s/%d.tmpf.%s", dir->s, file_idx, surfix);
  }
  s->l = strlen (s->s);

	return s->s;
}
