function diffmap_phase_retrieval_example1

close all;

opts.max_iterations=200;
example_opts.example_num=1;
example_opts.num_disks=50;
example_opts.N1=128;
example_opts.N2=128;
example_opts.oversamp=2;
random_seed=1;

rand_seed=clock; rand_seed=rand_seed(end);

rng(rand_seed);
opts.diffmap_method='AB';
example_opts.noise_level=0;
run_diffmap_phase_retrieval_example(example_opts,opts);
set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
pause(0.1);

rng(rand_seed);
opts.diffmap_method='BA';
example_opts.noise_level=0;
run_diffmap_phase_retrieval_example(example_opts,opts);
set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
pause(0.1);

rng(rand_seed);
opts.diffmap_method='AB';
example_opts.noise_level=0.05;
run_diffmap_phase_retrieval_example(example_opts,opts);
set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
pause(0.1);

rng(rand_seed);
opts.diffmap_method='BA';
example_opts.noise_level=0.05;
run_diffmap_phase_retrieval_example(example_opts,opts);
set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
pause(0.1);

end

