#include <sdf_read.h>
#include <string.h>

int sdf_read(SDF *g, const char *fname, int lev, double *time) {
	int j;
	int coord_count;
	int rc;

	if (lev < 0) {
		if(g->coords) free(g->coords);
		if(g->data) free(g->data);
		if(g->cnames[0]) free(g->cnames[0]);
		g->data=NULL;
		g->coords=NULL;
		g->cnames[0]=NULL;
		return(1);
	}
   strcpy(g->file,fname);

   rc = gft_read_rank(g->file,lev,&(g->dim));
	if (rc == 0) {return rc;}

   ivls(g->shape,1,max_dim);

   rc = gft_read_shape(g->file,lev,g->shape);	
	if (rc == 0) {return rc;}

	rc = gft_read_name(g->file,1,g->name);
	if (rc == 0) {return rc;}

	g->coord_size = 0;
	g->size = 1;

   for(j=0; j< g->dim; j++) {
		g->coord_size += g->shape[j];
		g->size *= g->shape[j];
	}	
	
	g->coords = vec_alloc(g->coord_size);
	g->data = vec_alloc(g->size);
	
	g->cnames = (char **)malloc(1*sizeof(char *));
	g->cnames[0] = (char *)malloc(6*sizeof(char));

	rc= gft_read_full(g->file,lev,g->shape,g->cnames[0],g->dim,time,
			g->coords, g->data);
	if (rc == 0) {return rc;}

	coord_count = 0;
	
	for(j=0;j<g->dim;j++) {
		g->bbox[j*2] = g->coords[coord_count];
      g->bbox[2*j+1] = g->coords[coord_count+g->shape[j]-1];
		coord_count += g->shape[j];
	}

   switch (g->dim) {
		case 3:
			g->x = g->coords;
			g->y = &(g->coords[g->shape[0]]);
			g->z = &(g->coords[g->shape[0]+g->shape[1]]);
			g->dx = g->x[1] - g->x[0];
			g->dy = g->y[1] - g->y[0];
			g->dz = g->z[1] - g->z[0];
			g->Nx = g->shape[0];
			g->Ny = g->shape[1];
			g->Nz = g->shape[2];
			g->xmin = g->bbox[0];
			g->xmax = g->bbox[1];
			g->ymin = g->bbox[2];
			g->ymax = g->bbox[3];
			g->zmin = g->bbox[4];
			g->zmax = g->bbox[5];
	
			break;
		case 2:
			g->x = g->coords;
			g->y = &(g->coords[g->shape[0]]);
			g->dx = g->x[1] - g->x[0];
			g->dy = g->y[1] - g->y[0];
			g->Nx = g->shape[0];
			g->Ny = g->shape[1];
			g->xmin = g->bbox[0];
			g->xmax = g->bbox[1];
			g->ymin = g->bbox[2];
			g->ymax = g->bbox[3];
			break;
		case 1:
			g->x = g->coords;
			g->dx = g->x[1] - g->x[0];
			g->Nx = g->shape[0];
			g->xmin = g->bbox[0];
			g->xmax = g->bbox[1];
			break;
		default:
			printf("sdf_read: **WARNING**, DIM not supported\n");
			break;
	}
  	
	return rc;
}
