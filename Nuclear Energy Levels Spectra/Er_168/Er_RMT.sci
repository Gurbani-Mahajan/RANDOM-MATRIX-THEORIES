clc
clear
//Uranium U 238

//loading and sorting energy levels
E_Er = [821.1685; ...
        1276.2716; ...
        1493.133; ...
        1848.354; ...
        1893.1; ...
        1930.391; ...
        1952.2; ...
        2193.19; ...
        2322.2; ...
        2450.5; ...
        2461.8; ...
        2538.1; ...
        2552.7; ...
        2741.9; ...
        2825; ...
        2878.9; ...
        2906; ...
        2933.44; ...
        2961.2; ...
        3009; ...
        3019.6; ...
        3042.3; ...
        3049.9; ...
        3055.95; ...
        3081.3; ...
        3098.4; ...
        3139.6; ...
        3172.5; ...
        3183.7; ...
        3194.4; ...
        3237.2; ...
        3269.4; ...
        3286.8; ...
        3342.9; ...
        3361.9; ...
        3429.2; ...
        3441.7; ...
        3451.6; ...
        3459.9; ...
        3471.6; ...
        3482.6; ...
        3493.3; ...
        3506.3; ...
        3515.7; ...
        3561.9; ...
        3617.6; ...
        3629.9; ...
        3720; ...
        3725.2; ...
        3740.4; ...
        3789.5; ...
        3808.5; ...
        3819.4; ...
        3861.9; ...
        3868.7; ...
        3876.3; ...
        3933; ...
        3964.9; ...
        4033.5; ...
        4055.9; ...
        4075.6];

sorted_E = gsort(E_Er, "g", "i");    // sorting energies ascending
sorted_E = sorted_E(:);  // ensuring that E_sorted is a column vector

//unfolding
function evals_unfold=unfold_spectrum(E_sorted)
    //removing duplicates
    E_unique=unique(E_sorted)
    n=length(E_unique)
    
    //finding N(E) (staircase func)
    N=1:n
    fit=polyfit(E_unique,N,6)
    evals_unfold=polyval(fit,E_unique)
    
    //plotting
    plot(E_unique,N,'o-b')
    plot(E_unique,evals_unfold,'r-')
endfunction

function spacing_u=find_nnsd(sorted_E)
    
sorted_E = sorted_E(:);
N = length(sorted_E);

//Compute average local spacing over sliding window
window = 10;  // Number of spacings to average over (adjustable)
s_local = zeros(N - 1, 1);
for i = 1:N-1
    low = max(1, i - floor(window/2));
    high = min(N - 1, i + floor(window/2))
    s_local(i) = mean(sorted_E(low+1:high+1) - sorted_E(low:high));
end

spacing_raw = diff(sorted_E);
spacing_u = spacing_raw ./ s_local;
spacing_u = spacing_u(find(spacing_u > 0 & ~isnan(spacing_u)));
spacing_u = spacing_u(:);
endfunction

//cdf (cumalative distribution function) of level spacings
function cdf=find_cdf(spacings)
    //assuming spacings are normalised P(s)=s
    spacings_sorted=gsort(spacings,'g','i') //sorting spacings
    n=length(spacings_sorted)
    cdf=(1:n)/n //since all spacings are normalised
endfunction

//sff (K)(spectral form factor) (K = |Z(T)|^2 where Z is partition function)
function sff=find_sff(e_sorted,T)
    evals_unfold=(e_sorted - min(e_sorted)) / (max(e_sorted) - min(e_sorted));
    n=length(evals_unfold)
    n_t=length(T)
    for i =1:n_t
    t=T(i)
    Z = sum(exp(2 * %pi * %i * evals_unfold * t))
    sff(i)= abs(Z)^2 / n;
    end
//sff=smooth(sff)
endfunction


//spectral rigidity Δ₃(L)
function Delta3_vals = compute_rigidity(E_sorted, L_vals)
    Delta3_vals = zeros(size(L_vals));
    N = length(E_sorted);

    for k = 1:length(L_vals)
        L = L_vals(k);
        sum_rigidity = 0;
        count = 0;

        for i = 1:N - L
            window_E = matrix(E_sorted(i:i+L-1), 1, L); //ensuring row vector
            N_x = 1:L;

            Xmat = [window_E' ones(L,1)];
            coeffs = Xmat \ N_x';
            N_fit = Xmat * coeffs;

            //corrected Δ₃ calculation: enforce non-negativity
            delta3_local = max(mean((N_x' - N_fit).^2),0); //0 if the local spectral rigidity is negative
            sum_rigidity = sum_rigidity + delta3_local;
            count = count + 1;
        end

        Delta3_vals(k) = sum_rigidity / count;
    end
endfunction

// finding chaoticity index
function CI = compute_chaos_index(L_vals, Delta3_emp)
    Delta3_poisson = matrix(L_vals / 15, 1, -1);
    Delta3_goe = matrix(log(L_vals) / (%pi^2), 1, -1);

    dist_to_poisson = mean(abs(Delta3_emp - Delta3_poisson));
    dist_to_goe = mean(abs(Delta3_emp - Delta3_goe));

    CI = dist_to_poisson / (dist_to_poisson + dist_to_goe);
endfunction

//finding empirical spectral observables
evals_unfold=unfold_spectrum(sorted_E)
spacing_u=find_nnsd(evals_unfold)
nnsd_sorted=gsort(spacing_u,'g','i')
cdf=find_cdf(spacing_u)
T=linspace(0,20,500)
sff=find_sff(sorted_E,T)
L_vals = floor(10.^linspace(log10(5), log10(20), 10));
Delta3_emp = compute_rigidity(sorted_E, L_vals);
CI = compute_chaos_index(L_vals, Delta3_emp);
disp("Chaoticity Index for Calcium: " + string(CI));

//theoretical values

//nnsd
s = 0:0.01:4
nnsd_goe = (%pi/2) * s .* exp(-(%pi/4) * s.^2);
nnsd_poisson = exp(-s); 

//cdf
y=linspace(0,5,300)
cdf_goe=1-exp(-(%pi*(y.^2))/4)
cdf_poisson=1-exp(-y)

//spectral rigidity
Delta3_poisson = L_vals / 15;
Delta3_goe = log(L_vals) / (%pi^2);


//comparing values by plotting


//nnsd
figure()
[counts,bins]=histplot(10,spacing_u)
//plot((bins(2:$-1))',P)
plot(s,nnsd_goe,'r-')
plot(s,nnsd_poisson, 'g-')
legend('Dy','GOE','Poisson')
xlabel('Nearest Neighbour Spacings')
ylabel('Frequency')
title('NNS histogram-Er 2+')

//cdf
figure()
plot(nnsd_sorted,cdf,'o-b')
plot(y,cdf_goe,'r-')
plot(y,cdf_poisson,'g-')
legend('Dy','GOE','Poisson')
xlabel('Normalised spacings')
ylabel('Cumalative Pro')
title('Level Spacings CDF- Er 2+')

//sff
figure()
plot(T,sff)
xlabel('T')
ylabel('K(t)')
title('Spectral Form Factor - Er 2+')

//spectral rigidity
figure()
plot(L_vals, (Delta3_emp), 'o-b') 
plot(L_vals, (Delta3_poisson), 'r--')
plot(L_vals, (Delta3_goe), 'g--')
legend("Dy", "Poisson", "GOE");
xlabel("Window Length L");
ylabel("log₁₀[Δ₃(L)]");
title("Spectral Rigidity- Er 2+");
