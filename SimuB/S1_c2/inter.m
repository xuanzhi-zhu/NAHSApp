%find the index of an event
indices=find((j-[0;j(1:end-1)])==1);

%corresponding t value
t_event=t(indices);

%inter-event intervals
t_inter=t_event-[0;t_event(1:end-1)];

figure(1)
plot(t_event(2:end),t_inter(2:end),'x')
set(gca, 'YScale', 'log')

%compute limit inferior inter-event intervals: PI=performance index
PI=zeros(size(t_event));
for i=1:1:numel(PI)
    PI(i)=1/i*sum(t_inter(1:i));
end

figure(2)
plot(t_event(2:end),PI(2:end),'x')
set(gca, 'YScale', 'log')

%0.721s