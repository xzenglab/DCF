function [s_cdf, s_prc] = save_evals(ind, funName, evals)

% save(sprintf('eval_%s_group_%d',funName, ind),'evals');
content = sprintf('eval for group %d: ',ind);


if ~isempty(evals{1,2})
    cdf_rates = evals{1,2};
    content = sprintf('%s %s',content,sprintf('cdf: %.3f',cdf_rates(100)));
    s_cdf = cdf_rates(100);
end

if and(~isempty(evals{2,2}),~isempty(evals{3,2}))
    pres = evals{3,2};
    recall = evals{2,2};
    prc = trapz(recall(1:100), pres(1:100));
    content = sprintf('%s %s',content,sprintf('AUPRC: %.3f',prc));
    s_prc = prc;
end

% subject = sprintf('%s evaluation got',funName);
% send_mail_upon_finished(subject, content,'18850544602@163.com');
disp('fun: save_evals done');
end

