function NewTimeSeries = SubtractClimatologyFromTimeSeries(TimeSeries)

time = length(TimeSeries);
NewTimeSeries = zeros(time,1);
for i=1:time
    Indices = mod(i,12):12:time;
    if mod(i,12)==0
        Indices = 12:12:time;
    end
    NewTimeSeries(i) = TimeSeries(i) - mean(TimeSeries(Indices));
end



end