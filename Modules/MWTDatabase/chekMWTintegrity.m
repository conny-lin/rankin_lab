function [result,pfiles] = chekMWTintegrity(pmwt)



%% check integrity of MWT files
[e,pset] =dircontent(pmwt,'*.set');
[b,pb] =dircontent(pmwt,'*.blobs');
if isempty(pb); [~,pb] =dircontent(pmwt,'*.blob'); end
[s,psum] =dircontent(pmwt,'*.summary');
[g,ppng] =dircontent(pmwt,'*.png');


if isempty(pset) || isempty(pb) || isempty(psum) || isempty(ppng)
   pfiles = dircontent(pmwt);
   result = 0;
else
    pfiles = [pset; pb; psum;ppng];
    result = 1;
end
    
end