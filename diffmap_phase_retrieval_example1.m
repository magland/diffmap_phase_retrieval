function diffmap_phase_retrieval_example1

%close all;

opts.max_iterations=5000;
example_opts.example_num=4;
example_opts.num_disks=50;
example_opts.N1=16;
example_opts.N2=16;
example_opts.oversamp=2;
example_opts.k=2;

rand_seed=5;
%rand_seed=clock; rand_seed=rand_seed(end);

if 1
rng(rand_seed);
opts.diffmap_method='AB';
opts.positivity=1;
opts.support=0;
opts.support2=0;
example_opts.noise_level=0;
run_diffmap_phase_retrieval_example(example_opts,opts);
set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
pause(0.1);
end;

if 0
rng(rand_seed);
opts.diffmap_method='BA';
opts.positivity=0;
opts.support=0;
opts.support2=1;
example_opts.noise_level=0;
run_diffmap_phase_retrieval_example(example_opts,opts);
set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
pause(0.1);
end;

% rng(rand_seed);
% opts.diffmap_method='AB';
% example_opts.noise_level=0;
% run_diffmap_phase_retrieval_example(example_opts,opts);
% set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
% pause(0.1);
% 
% rng(rand_seed);
% opts.diffmap_method='BA';
% example_opts.noise_level=0;
% run_diffmap_phase_retrieval_example(example_opts,opts);
% set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
% pause(0.1);
% 

% rng(rand_seed);
% opts.diffmap_method='AB';
% opts.eps_hack=0.01;
% example_opts.noise_level=0;
% run_diffmap_phase_retrieval_example(example_opts,opts);
% set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
% pause(0.1);
% 
% rng(rand_seed);
% opts.diffmap_method='BA';
% opts.eps_hack=0.1;
% example_opts.noise_level=0;
% run_diffmap_phase_retrieval_example(example_opts,opts);
% set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
% pause(0.1);

% 
% rng(rand_seed);
% opts.diffmap_method='AB';
% example_opts.noise_level=0.02;
% run_diffmap_phase_retrieval_example(example_opts,opts);
% set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
% pause(0.1);
% 
% rng(rand_seed);
% opts.diffmap_method='BA';
% example_opts.noise_level=0.02;
% run_diffmap_phase_retrieval_example(example_opts,opts);
% set(gcf,'name',sprintf('Method %s, Noise level %g',opts.diffmap_method,example_opts.noise_level));
% pause(0.1);

end

