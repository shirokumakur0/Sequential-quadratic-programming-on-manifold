function ImportGROPT
   GROPTBase=pwd;
   fprintf('Adding GROPT paths from %s...\n',GROPTBase);
   addpath(genpath(GROPTBase));
end
