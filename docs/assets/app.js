const candidates = {
  PACER0010: {
    generated: "PACERGEN01755", score: "0.7304", direction: -0.1251,
    note: "四条件分析已完成，但功能身份仍未解决，需要实验验证。"
  },
  PACER0027: {
    generated: "PACERGEN00123", score: "0.7307", direction: -0.0699,
    note: "四条件分析已完成；现有计算证据不足以判断 PAM 功能。"
  },
  PACER0073: {
    generated: "PACERGEN00094", score: "0.7041", direction: 0.1008,
    note: "正构激活核心出现重复同向变化（+0.603 / +0.211 / +0.261），仍属于待实验验证的动力学证据。"
  }
};

const nameEl = document.querySelector("#candidate-name");
const generatedEl = document.querySelector("#candidate-generated");
const scoreEl = document.querySelector("#candidate-score");
const directionEl = document.querySelector("#direction-value");
const directionBar = document.querySelector("#direction-bar");
const noteEl = document.querySelector("#candidate-note");

document.querySelectorAll(".candidate-tabs button").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll(".candidate-tabs button").forEach((item) => item.classList.remove("active"));
    button.classList.add("active");
    const id = button.dataset.id;
    const candidate = candidates[id];
    nameEl.textContent = id;
    generatedEl.textContent = candidate.generated;
    scoreEl.textContent = candidate.score;
    directionEl.textContent = candidate.direction >= 0
      ? "+" + candidate.direction.toFixed(4)
      : "−" + Math.abs(candidate.direction).toFixed(4);
    directionBar.style.left = (50 + candidate.direction * 100) + "%";
    noteEl.textContent = candidate.note;
  });
});
const observer = new IntersectionObserver((entries) => {
  entries.forEach((entry) => {
    if (entry.isIntersecting) {
      entry.target.classList.add("visible");
      observer.unobserve(entry.target);
    }
  });
}, { threshold: 0.12 });
document.querySelectorAll(".reveal").forEach((element) => observer.observe(element));

document.querySelector("#copy-command").addEventListener("click", async (event) => {
  const command = "cd submission\npython -m pip install -r requirements.txt\npython predict.py\npython verify_submission.py";
  try {
    await navigator.clipboard.writeText(command);
    event.currentTarget.textContent = "已复制";
  } catch {
    event.currentTarget.textContent = "请手动复制";
  }
});
