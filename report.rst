Поиск инфимума функционала
--------------------------

:Автор: Шачнев Дмитрий Алексеевич
:Группа: 411
:Дата: 2013-12-02

.. default-role:: math

Формализация задачи
===================

Требуется найти решение следующей задачи Лагранжа
при параметре `\alpha \in \{0,\;0.01,\;0.5,\;1.5,\;10.5\}`.

.. math::
   \left\{
   \eqalign{
    & \int\limits_0^{\pi / 2} u^2 dt + x_2^2(0) \to \inf \\
    & \dot x_1 = x_2 \\
    & \dot x_2 = u - \frac{x_1}{1 + \alpha t^2} \\
    & x_1(0) = 0 \\
    & x_1\left( \frac{\pi}{2} \right) = 1
   }
   \right.

Построение системы дифференциальных уравнений
=============================================

Функция Понтрягина:

.. math::
   H = p_1 x_2 + p_2 \left( u - \frac{x_1}{1 + \alpha t^2} \right) - \lambda_0 u^2

Уравнения Эйлера–Лагранжа:

.. math::
   \left\{
   \eqalign{
     \dot p_1 &= - \frac{\partial H}{\partial x_1} = \frac{p_2}{1 + \alpha t^2} \\
     \dot p_2 &= - \frac{\partial H}{\partial x_2} = -p_1
   }
   \right.

Функция Лагранжа:

.. math::
   \mathcal{L} = \int\limits_0^{\pi / 2} L dt + l

   L = \lambda_0 u^2 + p_1 (\dot x_1 - x_2) + p_2 \left( \dot x_2 - u + \frac{x_1}{1 + \alpha x_1^2} \right)
   \text{ — лагранжиан}

   l = \lambda_0 x_2^2(0) + \lambda_1 x_1(0) + \lambda_2 \left( x_1 \left( \frac{\pi}{2} \right) - 1 \right)
   \text{ — терминант}

Условия трансверсальности:

.. math::
   \left\{
   \eqalign{
    & p_1(0) = \frac{\partial\,l}{\partial x_1(0)} = \lambda_1 \\
    & p_1 \left( \frac{\pi}{2} \right) = \frac{\partial\,l}{\partial x_1 \left( \frac{\pi}{2} \right)}
      = -\lambda_2 \\
    & p_2(0) = \frac{\partial\,l}{\partial x_2(0)} = 2 \lambda_0 x_2(0) \\
    & p_2 \left( \frac{\pi}{2} \right) = \frac{\partial\,l}{\partial x_2 \left( \frac{\pi}{2} \right)} = 0
   }
   \right.

Условие стационарности:

.. math::
   2 \lambda_0 u - p_2 = 0

Если `\lambda_0 = 0`, то все множители Лагранжа равны одновременно нулю. Положим
`\lambda_0 = 1 / 2`. Тогда условие стационарности принимает вид `u = p_2`, и
итоговую систему можно записать в виде:

.. math::
   \left\{
   \eqalign{
    & \dot x_1 = x_2 \\
    & \dot x_2 = p_2 - \frac{x_1}{1 + \alpha t^2} \\
    & \dot p_1 = \frac{p_2}{1 + \alpha t^2} \\
    & \dot p_2 = -p_1 \\
    & x_1(0) = 0 \\
    & x_1\left( \frac{\pi}{2} \right) = 1 \\
    & p_2(0) = x_2(0) \\
    & p_2 \left( \frac{\pi}{2} \right) = 0
   }
   \right.

Решение системы при `\alpha = 0`
================================

При нулевом параметре `\alpha` система принимает вид:

.. math::
   \left\{
   \eqalign{
    & \dot x_1 = x_2 \\
    & \dot x_2 = p_2 - x_1 \\
    & \dot p_1 = p_2 \\
    & \dot p_2 = -p_1 \\
    & x_1(0) = 0 \\
    & x_1 \left( \frac{\pi}{2} \right) = 1 \\
    & p_2(0) = x_2(0) \\
    & p_2 \left( \frac{\pi}{2} \right) = 0
   }
   \right.

Эту систему можно решить аналитически. Из `\ddot p_1 = -p_1` имеем:

.. math::
    p_1(t) = C_1 \sin t + C_2 \cos t \\
    p_2(t) = C_1 \cos t - C_2 \sin t

Так как `p_2 \left( \frac{\pi}{2} \right) = 0`, то `C_2 = 0`, и

.. math::
    p_1(t) = C \sin t \\
    p_2(t) = C \cos t

Решение уравнения `\ddot x_1(t) = C \cos t - x_1(t)` даёт нам:

.. math::
   x_1(t) = A_1 \sin t + A_2 \cos t + \frac{C}{2} t \sin t

Из начальных условий `x_1(0) = 0,\; \dot x_1(0) = C,\; x_1 \left( \frac{\pi}{2} \right) = 1` находим:

.. math::
   A_2 = 0,\; A_1 = C,\; A_1 + \frac{C \pi}{4} = 1

   C = A_1 = \frac{4}{4 + \pi}

   \boxed{\displaystyle{ x_1(t) = \frac{4 + 2t}{4 + \pi} \sin t }}

Значение функционала при этом `x_1` равно:

.. math::
   \frac{4}{4 + \pi} \approx 0.56009915
