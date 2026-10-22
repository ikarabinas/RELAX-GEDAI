function test_RELAX_screen_robust_z()
    [z, c, s] = RELAX_screen_robust_z([1; 2; 3; 4; 100; NaN]);
    assert(c == 3 && abs(s - 1.4826) < 1e-12);
    assert(abs(z(5) - 97 / 1.4826) < 1e-9 && isnan(z(6)) && abs(z(1) + 2 / 1.4826) < 1e-9);

    [z, ~, s] = RELAX_screen_robust_z([5; 5; 5; 5; 9]);           % MAD = 0 -> mean-AD fallback
    assert(abs(s - 1.2533 * 0.8) < 1e-12 && abs(z(5) - 4 / (1.2533 * 0.8)) < 1e-9);

    [z, ~, s] = RELAX_screen_robust_z([2; 2; 2]);                  % no spread at all
    assert(s == 0 && all(isnan(z)));

    [z, c, s] = RELAX_screen_robust_z([NaN; NaN]);
    assert(all(isnan(z)) && isnan(c) && isnan(s));
    fprintf('test_RELAX_screen_robust_z: PASS\n');
end
