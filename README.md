AprilSAM: Real-time Smoothing and Mapping
===================================================
Copyright APRIL Lab.

https://april.eecs.umich.edu/

[AprilSAM: Real-time Smoothing and Mapping](https://april.eecs.umich.edu/papers/details.php?name=wang2018aprilsam)

[Video](https://april.eecs.umich.edu/public/users/xipengw/videos/AprilSAM.mp4)


INSTALL
=======

The default installation will place headers in /usr/local/include and
shared library in /usr/local/lib. It also installs a pkg-config script
into /usr/local/lib/pkgconfig.

    $ make
    $ sudo make install

USAGE
=====
### A basic AprilSAM application can be seen in `example/aprilsam_tutorial.c`.

Create a graph:

    april_graph_t *graph = april_graph_create();

Initialize optimize parameters:

    struct april_graph_cholesky_param *chol_param = calloc(1,sizeof(april_graph_cholesky_param_t));
    april_graph_cholesky_param_init(chol_param);
    chol_param->show_timing = 0; //
    chol_param->delta_xy = 0.1;
    chol_param->delta_theta = 0.1;
    chol_param->nthreshold = 100;

Add xyt node:

    double ninit[3]  = { 0,0,0 };
    double nstate[3] = { 0,0,0 };
    double ntruth[3] = { 0,0,0 };
    april_graph_node_t *node = april_graph_node_xyt_create(nstate, ninit, ntruth);
    node->UID = zarray_size(graph->nodes);
    zarray_add(graph->nodes, &node);

Add xyt factor:

    april_graph_node_t *na, *nb*;
    zarray_get(graph->nodes, aidx, &na);
    zarray_get(graph->nodes, bidx, &na);
    matd_t *W = matd_create(3,3);
    MATD_EL(W, 0, 0) = 1.0 / pow(0.1, 2);
    MATD_EL(W, 1, 1) = 1.0 / pow(0.1, 2);
    MATD_EL(W, 2, 2) = 1.0 / pow(to_radians(1), 2);
    april_graph_factor_t *factor = april_graph_factor_xyt_create(na->UID, nb->UID, z, NULL, W);
    zarray_add(graph->factors, &factor);
    matd_destroy(W);

Cholesky optimization:

    april_graph_cholesky(graph, chol_param);

AprilSAM increnmetal optimization:

    april_graph_cholesky_inc(graph, chol_param);

Graph Chi^2 error: (f(x)-z)^{T}\Sigma^{-1}(f(x)-z)

    april_graph_chi2(graph);

### An example of saving and loading graph file can be seen in `example/aprilsam_graph.c`

Save a graph:

    april_graph_save(april_graph_t *graph, const char *path)

Load a graph:

    april_graph_create_from_file(const char *path)

### An example of evaluting AprilSAM on Manhattan3500 dataset can be seen in `example/aprilsam_demo.c`

Test aprilsam optimization:

    ./aprilsam_demo --datapath ../data/M3500.txt

Test batch update only optimization:

    ./aprilsam_demo --datapath ../data/M3500.txt --batch_update_only
