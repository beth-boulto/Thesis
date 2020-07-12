%% ccode to calculate the expected number of mated parasites in a host with 
%M parasites in it for M going from 0 to 150 in increments of and plot the
%results

%set parasite vector
M = 0:150;
[~,n]=size(M);
%set mated female recording vector
F = zeros(size(M));
%loop over length of parasite vector
for i=1:n
   %if more than zero parasites in host 
   if M(i) > 0
       % initialise F1
       F1= 0;
       % for number of females from 1 to M-1 (need at least on efemale and
       % 1 male)
       for j=1:M(i)-1
           % add number of ways that number of females could be made from
           % that number of paparasites (M choose j), multiplied by the
           % number of females j
           F1=F1 + j*factorial(M(i))/(factorial(j)*factorial(M(i)-j));
       end
       % multiply by probability of number f females to get expected number
       % of females who will be mated - assumes as long as theres one male
       % all females will mate.
      F1= 2^(-M(i))*F1; 
       
   else
       %if no parasites - no females
       F1 = 0;
   end
    %record expected number of mated females
   F(i)=F1; 
       
end
Y = [M;F];
%plot expected number of females against total number of parasites
figure;
plot(M,F,'-.r','LineWidth', 4)
xlabel('Parasite population size, [M_i]')
ylabel('average number of mated females in population')