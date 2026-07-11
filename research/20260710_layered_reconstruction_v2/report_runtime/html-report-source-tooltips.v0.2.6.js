(() => {
  const RUNTIME_FLAG = "__dataAnalyticsSourceTooltipsInitialized";
  if (window[RUNTIME_FLAG]) return;
  window[RUNTIME_FLAG] = true;
  document.documentElement.setAttribute("data-source-tooltips-ready", "true");

  const TRIGGER_SELECTOR = ".source-tooltip";
  const CONTENT_SELECTOR = ".source-tooltip-content";
  const GAP = 8;
  const MOBILE_MAX_WIDTH = 600;
  const SCROLL_DISMISS_THRESHOLD = 24;
  const VIEWPORT_PADDING = 8;
  const FIGURE_TRIGGER_SELECTOR = ".source-figure > .source-tooltip";
  const FIGURE_BUTTON_ATTRIBUTE = "data-source-tooltip-figure-button";
  const FIGURE_STYLE_ID = "data-source-tooltip-figure-styles";

  let activeTrigger = null;
  let mobileTrayScrollY = 0;
  let placementFrame = 0;

  function viewportBounds() {
    const viewport = window.visualViewport;
    return {
      bottom: (viewport?.offsetTop ?? 0) + (viewport?.height ?? window.innerHeight),
      left: viewport?.offsetLeft ?? 0,
      right: (viewport?.offsetLeft ?? 0) + (viewport?.width ?? window.innerWidth),
      top: viewport?.offsetTop ?? 0,
    };
  }

  function usesMobileTray() {
    return (window.visualViewport?.width ?? window.innerWidth) <= MOBILE_MAX_WIDTH;
  }

  function upgradeFigureSourceButtons() {
    document.querySelectorAll(FIGURE_TRIGGER_SELECTOR).forEach((originalTrigger) => {
      if (!(originalTrigger instanceof HTMLElement) || originalTrigger instanceof HTMLButtonElement) return;
      const button = document.createElement("button");
      for (const attribute of Array.from(originalTrigger.attributes)) {
        if (attribute.name !== "tabindex" && attribute.name !== "role") {
          button.setAttribute(attribute.name, attribute.value);
        }
      }
      button.type = "button";
      button.setAttribute(FIGURE_BUTTON_ATTRIBUTE, "true");
      while (originalTrigger.firstChild) button.append(originalTrigger.firstChild);
      originalTrigger.replaceWith(button);
    });
  }

  function installFigureSourceStyles() {
    if (document.getElementById(FIGURE_STYLE_ID)) return;
    const style = document.createElement("style");
    style.id = FIGURE_STYLE_ID;
    style.textContent = `
      @media (max-width: 800px) {
        .chart-grid {
          grid-template-columns: minmax(0, 1fr);
          min-width: 0;
        }
        .source-figure {
          min-width: 0;
          max-width: 100%;
          overflow: hidden;
        }
        .source-figure .chart-wrap {
          max-width: 100%;
          overflow-x: auto;
        }
      }
    `;
    document.head.append(style);
  }

  function clampAxis(value, size, start, end) {
    const minimum = start + VIEWPORT_PADDING;
    const maximum = Math.max(minimum, end - size - VIEWPORT_PADDING);
    return Math.min(Math.max(value, minimum), maximum);
  }

  function placeTooltip(trigger) {
    if (usesMobileTray() || !(trigger instanceof HTMLElement)) return;
    const tooltip = trigger.querySelector(CONTENT_SELECTOR);
    if (!(tooltip instanceof HTMLElement)) return;

    const triggerRect = trigger.getBoundingClientRect();
    const tooltipRect = tooltip.getBoundingClientRect();
    const viewport = viewportBounds();
    const isFigureTooltip = trigger.parentElement?.classList.contains("source-figure");
    const centeredLeft = triggerRect.left + (triggerRect.width - tooltipRect.width) / 2;
    const aboveTop = triggerRect.top - tooltipRect.height - GAP;
    const belowTop = triggerRect.bottom + GAP;
    const preferredTop = isFigureTooltip || aboveTop < viewport.top + VIEWPORT_PADDING
      ? belowTop
      : aboveTop;
    const left = clampAxis(centeredLeft, tooltipRect.width, viewport.left, viewport.right);
    const top = clampAxis(preferredTop, tooltipRect.height, viewport.top, viewport.bottom);

    tooltip.style.setProperty("--source-tooltip-left", `${Math.round(left)}px`);
    tooltip.style.setProperty("--source-tooltip-top", `${Math.round(top)}px`);
  }

  function triggerValue(trigger, tooltip) {
    return Array.from(trigger.childNodes)
      .filter((node) => node !== tooltip)
      .map((node) => node.textContent ?? "")
      .join("")
      .trim();
  }

  function prepareMobileTray(trigger) {
    const tooltip = trigger.querySelector(CONTENT_SELECTOR);
    if (!(tooltip instanceof HTMLElement)) return false;

    let heading = tooltip.querySelector(".source-tooltip-heading");
    if (!(heading instanceof HTMLElement)) {
      heading = document.createElement("span");
      heading.className = "source-tooltip-heading";
      tooltip.prepend(heading);
    }
    const value = triggerValue(trigger, tooltip);
    heading.textContent = value ? `Source for ${value}` : "Source details";
    return true;
  }

  function closeMobileTray({ restoreFocus = false } = {}) {
    const trigger = activeTrigger;
    if (!(trigger instanceof HTMLElement)) return;
    trigger.removeAttribute("data-source-tooltip-open");
    trigger.setAttribute("aria-expanded", "false");
    activeTrigger = null;
    if (restoreFocus) trigger.focus({ preventScroll: true });
  }

  function openMobileTray(trigger) {
    if (!(trigger instanceof HTMLElement) || !prepareMobileTray(trigger)) return;
    if (activeTrigger && activeTrigger !== trigger) closeMobileTray();
    activeTrigger = trigger;
    mobileTrayScrollY = window.scrollY;
    trigger.setAttribute("data-source-tooltip-open", "true");
    trigger.setAttribute("aria-expanded", "true");
  }

  function toggleMobileTray(trigger) {
    if (activeTrigger === trigger && trigger.getAttribute("data-source-tooltip-open") === "true") {
      closeMobileTray();
    } else {
      openMobileTray(trigger);
    }
  }

  function syncSemantics() {
    const isMobile = usesMobileTray();
    document.querySelectorAll(TRIGGER_SELECTOR).forEach((trigger) => {
      if (!(trigger instanceof HTMLElement)) return;
      if (isMobile) {
        if (!trigger.hasAttribute("role")) {
          trigger.setAttribute("role", "button");
          trigger.setAttribute("data-source-tooltip-runtime-role", "true");
        }
        trigger.setAttribute("aria-expanded", trigger.getAttribute("data-source-tooltip-open") === "true" ? "true" : "false");
      } else {
        if (trigger.getAttribute("data-source-tooltip-runtime-role") === "true") {
          trigger.removeAttribute("role");
          trigger.removeAttribute("data-source-tooltip-runtime-role");
        }
        trigger.removeAttribute("aria-expanded");
      }
    });
  }

  function schedulePlacement(trigger = activeTrigger) {
    if (!(trigger instanceof HTMLElement)) return;
    activeTrigger = trigger;
    if (placementFrame) window.cancelAnimationFrame(placementFrame);
    placementFrame = window.requestAnimationFrame(() => {
      placementFrame = 0;
      placeTooltip(trigger);
    });
  }

  function eventTrigger(event) {
    const target = event.target;
    return target instanceof HTMLElement ? target.closest(TRIGGER_SELECTOR) : null;
  }

  function handleClick(event) {
    if (!usesMobileTray()) return;
    const target = event.target;
    if (!(target instanceof HTMLElement)) return;
    if (target.closest(CONTENT_SELECTOR)) return;
    const trigger = target.closest(TRIGGER_SELECTOR);
    if (trigger instanceof HTMLElement) {
      event.preventDefault();
      event.stopPropagation();
      toggleMobileTray(trigger);
      return;
    }
    closeMobileTray();
  }

  function handleKeydown(event) {
    if (!usesMobileTray()) return;
    if (event.key === "Escape") {
      closeMobileTray({ restoreFocus: true });
      return;
    }
    if (event.key !== "Enter" && event.key !== " ") return;
    const trigger = eventTrigger(event);
    if (!(trigger instanceof HTMLElement)) return;
    event.preventDefault();
    toggleMobileTray(trigger);
  }

  function handleScroll() {
    if (usesMobileTray()) {
      if (activeTrigger && Math.abs(window.scrollY - mobileTrayScrollY) >= SCROLL_DISMISS_THRESHOLD) {
        closeMobileTray();
      }
      return;
    }
    schedulePlacement();
  }

  function handleViewportChange() {
    syncSemantics();
    if (usesMobileTray()) {
      if (activeTrigger?.getAttribute("data-source-tooltip-open") !== "true") activeTrigger = null;
      return;
    }
    if (activeTrigger?.getAttribute("data-source-tooltip-open") === "true") {
      closeMobileTray();
      return;
    }
    schedulePlacement();
  }

  function setup() {
    upgradeFigureSourceButtons();
    installFigureSourceStyles();
    syncSemantics();
    document.addEventListener("pointerover", (event) => {
      if (!usesMobileTray()) schedulePlacement(eventTrigger(event));
    }, true);
    document.addEventListener("focusin", (event) => {
      if (!usesMobileTray()) schedulePlacement(eventTrigger(event));
    });
    document.addEventListener("click", handleClick, true);
    document.addEventListener("keydown", handleKeydown);
    document.addEventListener("scroll", handleScroll, true);
    window.addEventListener?.("resize", handleViewportChange);
    window.visualViewport?.addEventListener("resize", handleViewportChange);
    window.visualViewport?.addEventListener("scroll", handleScroll);

    const currentTrigger = document.querySelector(`${TRIGGER_SELECTOR}:hover, ${TRIGGER_SELECTOR}:focus`);
    if (!usesMobileTray() && currentTrigger instanceof HTMLElement) schedulePlacement(currentTrigger);
  }

  setup();
})();
