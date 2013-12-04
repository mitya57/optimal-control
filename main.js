 /***********************************\
/   Finding infimum of a functional   \
\    Author: Dmitry Shachnev, 2013    /
 \***********************************/

//
//  These are constants you may want to modify according to
//  your needs.
//

var EPS = 1e-12;
var TAU = 2e-5;

var defaults = {
    X1: function(t) {
        return (4 + 2 * t) * Math.sin(t) / (4 + Math.PI);
    },
    X2: function(t) {
        return (2 * Math.sin(t) + 4 * Math.cos(t) +
                2 * t * Math.cos(t)) / (4 + Math.PI);
    },
    P1: function(t) {
        return 4 * Math.sin(t) / (4 + Math.PI);
    },
    P2: function(t) {
        return 4 * Math.cos(t) / (4 + Math.PI);
    }
};

var derivatives = {
    X1: function(state, t, alpha) {
        return state.X2;
    },
    X2: function(state, t, alpha) {
        return state.P2 - state.X1 / (1 + alpha * t * t);
    },
    P1: function(state, t, alpha) {
        return state.P2 / (1 + alpha * t * t);
    },
    P2: function(state, t, alpha) {
        return -state.P1;
    }
};

var alphas = [0, 0.01, 0.5, 1.5, 10.5];

//
//  This code is mostly task-independent.
//

var integrate = function(func, a, b, tau) {
    var result = 0;
    for (var t = a; t + tau < b; t += tau) {
        result += (func(t) + func(t + tau)) * tau / 2;
    }
    return result;
};

var square = function(func) {
    return function(t) {
        return func(t) * func(t);
    };
};

var runge_kutta_diff = function(state, t, der, alpha, tau) {
    /* FIXME: actually implement Runge-Kutta here */
    return der(state, t, alpha) * tau;
};

var norm = function(v) {
    return Math.sqrt(v[0] * v[0] + v[1] * v[1]);
};

var get_values = function(P1_0, P2_0, alpha, tau, point, cache) {
    var state = {
        X1: 0,
        X2: P2_0,
        P1: P1_0,
        P2: P2_0
    };
    var states = [];
    for (var t = 0; t < point; t += tau) {
        state.X1 += runge_kutta_diff(state, t, derivatives.X1, alpha, tau);
        state.X2 += runge_kutta_diff(state, t, derivatives.X2, alpha, tau);
        state.P1 += runge_kutta_diff(state, t, derivatives.P1, alpha, tau);
        state.P2 += runge_kutta_diff(state, t, derivatives.P2, alpha, tau);
        if (cache) {
            states.push({
              X1: state.X1,
              X2: state.X2,
              P1: state.P1,
              P2: state.P2
            });
        }
    }
    return cache ? states : state;
};

var get_boundary_diff = function(P1_0, P2_0, alpha, tau) {
    var values = get_values(P1_0, P2_0, alpha, tau, Math.PI / 2, false);
    return [values.X1 - 1, values.P2];
};

var reverse_matrix = function(matrix) {
    var det = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    return [
       matrix[3] / det,
      -matrix[1] / det,
      -matrix[2] / det,
       matrix[0] / det
    ];
};

var function_from_cache = function(cache, fname, tau) {
    return function(t) {
        return cache[Math.floor(t / tau)][fname];
    }
};

var find_infimum = function(alpha, eps, tau) {
    if (!alpha) {
        return defaults;
    }
    var P1_0 = 0;
    var P2_0 = 4 / (4 + Math.PI);
    var diff0, diff1, diff2, pd, rd, iter = 0;
    var P1_0_new, P2_0_new, gamma, diff_norm, diff_norm_prev;
    var delta = eps;
    do {
        diff0 = get_boundary_diff(P1_0, P2_0, alpha, tau);
        diff1 = get_boundary_diff(P1_0 + delta, P2_0, alpha, tau);
        diff2 = get_boundary_diff(P1_0, P2_0 + delta, alpha, tau);
        pd = [
          (diff1[0] - diff0[0]) / delta,
          (diff2[0] - diff0[0]) / delta,
          (diff1[1] - diff0[1]) / delta,
          (diff2[1] - diff0[1]) / delta
        ];
        rd = reverse_matrix(pd);
        gamma = 1;
        if (iter) {
            diff_norm_prev = diff_norm;
            do {
                P1_0_new = P1_0 - gamma * (rd[0] * diff0[0] + rd[1] * diff0[1]);
                P2_0_new = P2_0 - gamma * (rd[2] * diff0[0] + rd[3] * diff0[1]);
                diff1 = get_boundary_diff(P1_0_new, P2_0_new, alpha, tau);
                /* Fedorenko norm */
                diff_norm = Math.sqrt(
                  diff1[0] * diff1[0] / (pd[0] * pd[0] + pd[1] * pd[1]) +
                  diff1[1] * diff1[1] / (pd[2] * pd[2] + pd[3] * pd[3])
                );
                gamma /= 2;
            } while (diff_norm > diff_norm_prev);
            P1_0 = P1_0_new;
            P2_0 = P2_0_new;
        } else {
            diff_norm = norm(diff0);
        }
        //console.log("norm of diff0 is " + norm(diff0));
        ++iter;
    } while (diff_norm > eps);
    console.log('found optimal boundary conditions: ' + P1_0 + ', ' + P2_0);
    var cache = get_values(P1_0, P2_0, alpha, tau, Math.PI / 2, true);
    return {
        X1: function_from_cache(cache, 'X1', tau),
        X2: function_from_cache(cache, 'X2', tau),
        P2: function_from_cache(cache, 'P2', tau)
    };
};

var process = function(x1, x2, u, alpha, tau) {
    var integral = integrate(square(u), 0, Math.PI / 2, tau);
    //console.log("integral = " + integral);
    var terminal = x2(0) * x2(0);
    //console.log("terminal = " + terminal);
    return integral + terminal;
};

var print_result = function(alpha, eps, tau) {
    console.log('***** alpha = ' + alpha + ' *****');
    var funcs = find_infimum(alpha, eps, tau);
    var result = process(funcs.X1, funcs.X2, funcs.P2, alpha, tau);
    console.log('result = ' + result);
};

alphas.forEach(function(alpha) {
    print_result(alpha, EPS, TAU);
});
