%function alpha = line_search(func, d, x0, f_min, alpha1, c1, c2)
%%%%%%%%%for testing%%%%%%%%%
function alpha = line_search(d, x0, f_min, alpha1, c1, c2)
func = @(x)func_trial(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_plot = true;
%LINE_SEARCH conduct line search under strong wolfe condition for a
%multi-variable function on a certain direction
%input: 
% func: the multi-variable function, return value and gradient
% d: searching direction
% x0: starting point
% f_min: estimation of minimum value of func
% alpha1: initial value of alpha
% c1, c2: coefficients for wolfe conditions, 0<c1<c2<1, typically c1~0,
% c2~1 for loose search.
%output:
% alpha: the step length satisfying strong wolfe condition:
% (1)  func(x0+alpha*d)<=func(x0)+c1*alpha*d'*grad_func(x0+alpha*d),
% (2)  d'*grad_func(x0+alpha*d)<=-|c2*d'*grad_func(x0)|
% see pp. 59-61 of Nocedal, J. & Wright, S. J., 1999. Numerical optimization.
alpha_pre = 0;
[f0, g0] = func(x0);
alpha_max = - (f0 - f_min) / (c1 * d' * g0);
alpha_cur = alpha_max / 1000;
if (alpha1 < alpha_cur)
    alpha_cur = alpha1;
else
    [f1, g1] = func(x0 + alpha1 * d);
    if ((f1 > f0 + c1 * alpha1 * (d' * g0)) && (abs(d' * g1) <= - c2 * (d' * g0)))
        return
    end
end
f_pre = f0;
if show_plot
    figure(1);
    alphas = 0: alpha_max / 1000: alpha_max;
    [fs, ~] = func(bsxfun(@plus, x0, alphas * d));
    lines = bsxfun(@plus, f0, c1 * alphas * (d' * g0));
    plot(alphas, fs);
    hold on
    plot(alphas, lines);
end
while (alpha_pre < alpha_max)
    if (alpha_cur > alpha_max)
        alpha_cur = alpha_max;
    end
    [f_cur, g_cur] = func(x0 + alpha_cur * d);
    if show_plot
        figure(1);
        hold on
        plot(alpha_cur, f_cur, 'o');
    end
    if ((f_cur > f0 + c1 * alpha_cur * (d' * g0)) || (f_cur >= f_pre))
        alpha = zoom(func, alpha_pre, alpha_cur, d, x0, c1, c2, f0, g0);
        return
    end
    if (abs(d' * g_cur) <= - c2 * (d' * g0))
        alpha = alpha_cur;
        return
    end
    if (d' * g_cur >= 0)
        alpha = zoom(func, alpha_cur, alpha_pre, d, x0, c1, c2, f0, g0);
        return
    end
    alpha_pre = alpha_cur;
    alpha_cur = 200 * alpha_pre;
    f_pre = f_cur;
end
alpha = alpha_max;
end

function alpha = zoom(func, alpha_lo, alpha_hi, d, x0, c1, c2, f0, g0)
show_plot = true;
[f_lo, g_lo] = func(x0 + alpha_lo * d);
[f_hi, g_hi] = func(x0 + alpha_hi * d);
max_iter = 20;
N_iter = 0;
error_tol = 1e-5;
seg_len = alpha_hi - alpha_lo;
if show_plot
    figure(2);
    alphas = alpha_lo: (alpha_hi - alpha_lo) / 1000: alpha_hi;
    [fs, ~] = func(bsxfun(@plus, x0, alphas * d));
    lines = bsxfun(@plus, f0, c1 * alphas * (d' * g0));
    plot(alphas, fs);
    hold on
    plot(alphas, lines);
end
while ((abs(seg_len) > error_tol) && (N_iter <= max_iter))
    eta = 3 * (f_hi - f_lo) - d' * (2 * g_lo + g_hi) * seg_len;
    xi = d' * (g_lo + g_hi) * seg_len - 2 * (f_hi - f_lo);
    if (abs(xi) < error_tol)
        alpha = alpha_lo + 0.5 * (d' * g_lo) * seg_len * seg_len / ((d' * g_lo) * seg_len + f_hi - f_lo);
    else
        Delta = sqrt(eta * eta - 3 * xi * (d' * g_lo) * seg_len);
        z1 = (-eta + Delta) / (3 * xi);
        z2 = (-eta - Delta) / (3 * xi);
        if (min(z1, z2) < 0)
            z = max(z1, z2);
        else
            z = min(z1, z2);
        end
        alpha = alpha_lo + z * seg_len;
    end
    [f_here, g_here] = func(x0 + alpha * d);
    if show_plot
        figure(2);
        hold on
        plot(alpha, f_here, 'o');
    end
    if ((f_here > f0 + c1 * alpha * (d' * g0)) || (f_here >= f_lo))
        %alpha_hi = alpha;
        f_hi = f_here;
        g_hi = g_here;
        seg_len = alpha - alpha_lo;
    else
        if (abs(d' * g_here) <= - c2 * (d' * g0))
            return
        end
        if (d' * g_here * seg_len >= 0)
            alpha_hi = alpha_lo;
            f_hi = f_lo;
            g_hi = g_lo;
        end
        f_lo = f_here;
        g_lo = g_here;
        alpha_lo = alpha;
        seg_len = alpha_hi - alpha_lo;           
    end
end
end

function [f,g] = func_trial(x)
f = cos(x);
g = -sin(x);
end

